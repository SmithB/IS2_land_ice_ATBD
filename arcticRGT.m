function [RGT, laser] = arcticRGT(roi, sample_rate)

% ROI bounds
roi_xlim = [min([roi.X]) max([roi.X])];
roi_ylim = [min([roi.Y]) max([roi.Y])];

% extract satellite ground track
foo = shaperead(fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data','RGTs', 'arctic_rgt.shp'));
n = length(foo);
[RGT(1:n).X] = deal([]);
[RGT(1:n).Y] = deal([]);
[RGT(1:n).x_RGT] = deal([]);
[RGT(1:n).trackID] = deal([]);
%% PARFOR
for j = 1:n
    RGT(j).X = foo(j).X;
    RGT(j).Y = foo(j).Y;
    RGT(j).x_RGT = ...
        cumsum(sqrt(([0 RGT(j).X(2:end) - ...
        RGT(j).X(1:end-1)]).^2 + ...
        ([0 RGT(j).Y(2:end) - ...
        RGT(j).Y(1:end-1)]).^2));
    RGT(j).trackID = str2double(foo(j).Name);
end

% extract reference ground tracks
[laser0(1:6).track] = deal([]);

% for lasers 1-6
%!!simplifying approximation!! - set accross track ditance to RGT as fixed
%xt = [3290 3380 -45 45 -3380 3290];
xt = [-3290 -3380; 45 -45; 3380 3290];

%% PARFOR
for i = 1:6
    foo = shaperead(fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data','RGTs', ['arctic_rgt_laser' num2str(i) '.shp']));
  
    %!!simplifying approximation!! - set track distance along track
    %starting arbitrarily from start of line  segement
    for j = 1:length(foo)
        laser0(i).track(j).x_RGT = ...
            cumsum(sqrt(([0 foo(j).X(2:end) - ...
            foo(j).X(1:end-1)]).^2 + ...
            ([0 foo(j).Y(2:end) - ...
            foo(j).Y(1:end-1)]).^2));
        laser0(i).track(j).y_RGT = ...
            zeros(size(laser0(i).track(j).x_RGT)) + xt(i);
        laser0(i).track(j).trackID = str2double(foo(j).Name);
        laser0(i).track(j).X = foo(j).X;
        laser0(i).track(j).Y = foo(j).Y;
    end 
end

% find sections of RGT that cross ROI
% and delete all data that does not intersect buffered ROI
buffer = 500E3;
delIdx = false(size(RGT,1));

figure; plot(roi_xlim([1 1 2 2 1]), roi_ylim([1 2 2 1 1]),'r'); hold on
%% PARFOR
for j = 1:length(RGT)
    S = RGT(j);
     
    in_roi =  S.X >= (roi_xlim(1) - buffer) & S.X <= (roi_xlim(2) + buffer) & ...
        S.Y >= (roi_ylim(1) - buffer) &  S.Y <= (roi_ylim(2) + buffer);
    delIdx(j) = sum(in_roi) <= 1;
    
    if ~delIdx(j)
        % only include segments that intersect region
        S = removeStucData(S,~in_roi);
        
        
        % interpolate along track at 10 m spacing
        d_dense = S.x_RGT(1):10:S.x_RGT(end);
        
        S = interpStructData(S, S.x_RGT, d_dense);
        
        % trim to roi box
        in_roi =  S.X >= (roi_xlim(1)) & S.X <= (roi_xlim(2)) & ...
            S.Y >= (roi_ylim(1)) &  S.Y <= (roi_ylim(2));
        
        % only include segments that intersect region
        S = removeStucData(S,~in_roi);
        
        % trim to roi polygon
        in_roi = inpolygon(S.X, S.Y, roi.X(~isnan(roi.X)), roi.Y(~isnan(roi.Y)));
        
        delIdx(j) = sum(in_roi) <= 1;
        if ~delIdx(j)
            % only include segments that intersect region
            S = removeStucData(S,~in_roi);
            
            % sample at rate of sensor
            d_dense = S.x_RGT(1):sample_rate:S.x_RGT(end);
            RGT(j) = interpStructData(S, S.x_RGT,d_dense);
        end
        plot(RGT(j).X, RGT(j).Y,'kx');
    end
end
m0 = length(RGT);
RGT(delIdx) = [];
trackID = [RGT(:).trackID];

% find matching track ID's
for i = 1:6
    m = length(laser0(i).track);
    if m == m0
        laser(i).track = laser0(i).track(~delIdx);
    else
        % not all RGT have the same number of traks... not sure why
        idx = ismember([laser0(i).track(:).trackID],trackID);
        laser(i).track = laser0(i).track(idx);
    end
end

% find sections of RGT that cross ROI
% and delete all data that does not intersect buffered ROI

for i = 1:6
    delIdx = false(length(laser(i).track),1);
    for j = 1:length(laser(i).track)
        S = laser(i).track(j);
        
        buffer = 1000E3;
        in_roi =  S.X >= (roi_xlim(1) - buffer) & S.X <= (roi_xlim(2) + buffer) & ...
            S.Y >= (roi_ylim(1) - buffer) &  S.Y <= (roi_ylim(2) + buffer);
        
        % only include segments that intersect region
        S = removeStucData(S,~in_roi);
        if isempty(S.x_RGT) || length(S.x_RGT)<=1
            delIdx(j) = true;
            continue
        end
        
        % interpolate along track at 10 m spacing
        d_dense = S.x_RGT(1):10:S.x_RGT(end);
        
        S = interpStructData(S, S.x_RGT, d_dense);
        
        buffer = 100E3;
        % trim to roi box
        in_roi =  S.X >= (roi_xlim(1) - buffer) & S.X <= (roi_xlim(2) + buffer) & ...
            S.Y >= (roi_ylim(1) - buffer) &  S.Y <= (roi_ylim(2) + buffer);
         
        % only include segments that intersect region
        S = removeStucData(S,~in_roi);
        if  isempty(S.x_RGT) || length(S.x_RGT)<=1
            delIdx(j) = true;
            continue
        end
        
        % sample at rate of sensor
        d_dense = S.x_RGT(1):sample_rate:S.x_RGT(end);
        S = interpStructData(S, S.x_RGT,d_dense);
        
        % find closest point to satellite ground track
        
        % closest start and end point to RGT to make lines parallel
        idxRGT = S.trackID == trackID;
        ds = sqrt((S.X - RGT(idxRGT).X(1)).^2 + (S.Y - RGT(idxRGT).Y(1)).^2);
        de = sqrt((S.X - RGT(idxRGT).X(end)).^2 + (S.Y - RGT(idxRGT).Y(end)).^2);
        ds = find(min(ds) == ds);
        de = find(min(de) == de);
        
        sIdx = min(ds,de);
        eIdx = max(ds,de);
        
%         plot(RGT(j).X,RGT(j).Y, '.g'); hold on
%         plot(S.X,S.Y, '.r')
%         plot(S.X(sIdx),S.Y(sIdx), '.k', 'MarkerSize', 5)
%         plot(S.X(eIdx),S.Y(eIdx), '.k', 'MarkerSize', 5)

        % make sure the satellite and laser RGTs are the same length
        d = length(RGT(idxRGT).X) - (eIdx-sIdx+1);
        de = ceil(d/2);
        ds = d-de;
        sIdx = sIdx-ds;
        eIdx = eIdx+de;
        parallelIdx = false(size(S.X));
        parallelIdx(sIdx:eIdx) = true;
        laser(i).track(j) = removeStucData(S,~parallelIdx);
        laser(i).track(j).x_RGT = RGT(idxRGT).x_RGT; 
    end
    laser(i).track(delIdx) = [];
end






