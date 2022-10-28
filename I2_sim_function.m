function I2_sim_function(tileObjectID,acquisition_step)

tic
time0 = now;

addpath('/Volumes/ice2/nick/34_IceSat2/SyntheticCode_Alex')

%% User options
% plot intermediate output
plot_flag = false;
simulated_topography = 'checker'; % required if run_type = 3 ['peaks', 'checker']
single_track = false; % only simulates the longest track in the roi
parallel_flag = false;
ben_flag = 1; % If on Bens Computer

% specify ROI shapefile name
% ###### must fall entirely within a single ArcticDEM tile footprint ######
roi_shapefile = 'ROI_NEGIS.shp'; % required if run_type = 1
% #########################################################################



% %%%%%%%% TileID's for Jakobshavn
tileObjectID_opts = [15327;15328;14135;14133;7220;7217;7153;7151;15329;15326;14132;14134;7218;7219;7154;7152; ...
    13621;13620;15028;15030;7221;7222;7122;7121;13618;13619;15029;15027;7223;7224;7120;7119; ...
    14522;14520;11998;11997;7172;7173;7140;7141;14523;14521;12000;11999;7174;7175;7139;7142; ...
    12869;12867;13782;13781;7159;7161;7131;7134;12868;12866;13784;13783;7160;7162;7132;7133];


load /Volumes/ice2/nick/34_IceSat2/SyntheticCode_Alex/SignalAdditions.mat
load /Volumes/ice2/nick/34_IceSat2/SyntheticCode_Alex/Mispointing.mat

kk = find(tileObjectID_opts == tileObjectID);


%%%%%%%%%%%%%%%%%%%% Set the arcticDEM directory

ad_local_dir = '/Volumes/ice2/nick/34_IceSat2/SyntheticCode_Alex/Data/ArcticDEM/tiles/';

if exist(fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', ['atl03_sim2_',sprintf('%0.2d',acquisition_step)])) == 0
    mkdir(fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', ['atl03_sim2_',sprintf('%0.2d',acquisition_step)]));
end

if exist(fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', ['atl06_sim2_',sprintf('%0.2d',acquisition_step)])) == 0
    mkdir(fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', ['atl06_sim2_',sprintf('%0.2d',acquisition_step)]));
end




%% simulation parameters

% params: a sturct with simulation parameter fields:
params.roughness = .3; % surface roughness within footprint [m]
params.sample_rate = 0.7; % along track [m]
params.pulse_frequency = 10E3; % Hz [s^-1]
params.ground_speed = params.sample_rate * params.pulse_frequency;
params.sigma_x = 7.2; % footprint width [m]
params.sigma_pulse = 1.6e-9;
params.H_window = 600;% size of surface window [m]
params.c = 3e8; % speed of light, m/s
params.t_dead = 3.2e-9; % 3.2e-9 dead time, in seconds
params.pulses_per_seg = 57;

% beam specific parameters
params.NoiseRate = 8e6;% noise rate [Hz]
params.N_channels = 16; % number of channels in the detector
params.N_per_pulse = 12; % mean number of detected photons in the return pulse

% now for weak beam
params(2) = params;
params(2).N_channels = 4;
params(2).NoiseRate = params(2).N_channels/ params(1).N_channels * params(1).NoiseRate;
params(2).N_per_pulse = params(2).N_channels/ params(1).N_channels * params(1).N_per_pulse;

% file paths

wf_file = fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', 'Parameters', 'WF_est.mat');
snr_f_file = fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', 'Tables', 'SNR_F_table.h5');
ad_tile_idx_file = fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', 'ArcticDEM', 'ArcticDEM_Tile_Index_2017sept06.shp');

strongBeam = [1,3,5];
weakBeam = [2,4,6];

% struct giving transmit-pulse power as a function of time, with fields t and p
A3 = load(wf_file);
params(1).WF = A3.WF;
params(2).WF = A3.WF;


%% deteremine runName

runName = num2str(tileObjectID);


%% Define ROI with sapefile & find ArcticDEM tile
% open up a parpool for parallel processing [optional]
% if license('test', 'image_toolbox') == 1
%     foo = gcp;
%     if parallel_flag && isempty(foo)
%         parpool(3)
%     end
% end

% update shapefile reader... this only needs to be run the first time the
% program is run... but no harm running it very time
% shapereadUpdate
addpath(pwd)

ad_tile = shaperead(ad_tile_idx_file);
in_roi = tileObjectID == [ad_tile(:).objectid];
ad_tile = ad_tile(in_roi);

roi = ad_tile;
roi_xlim = [min([roi.X]) max([roi.X])];
roi_ylim = [min([roi.Y]) max([roi.Y])];

if length(ad_tile) > 1
    error('ROI crosses multiple ArctciDEM tile footprints which is not allowed, please reposition ROI')
elseif isempty(ad_tile)
    %error('ROI does not intersect an ArctciDEM tile footprint, please reposition ROI')
    dem_flag = 1;
else
    dem_flag = 0;
end

%% Download ArcticDEM tiles [if they do not exist already]
if dem_flag == 0
    for i = 1:length(ad_tile)
        fileurl = ad_tile(i).fileurl;
        [~,fName,ext] = fileparts(fileurl);
        
        % check if folder already exists
        if exist(fullfile(ad_local_dir, fName), 'dir') ~= 7
            
            if ben_flag == 1
                ad_local_dir = fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', 'ArcticDEM', 'tiles');
                fileurl = ad_tile(i).fileurl;
                [~,fName,ext] = fileparts(fileurl);
            end
            
            % make directory if it doesn't exist
            if exist(ad_local_dir, 'dir') ~= 7
                mkdir(ad_local_dir)
            end
            fprintf('Downloading ArcticDEM %s: %.1fmin\n', runName, (now-time0)*24*60)
            
            local_path_tar = fullfile(ad_local_dir, [fName,ext]);
            
            % save in tempfile then rename
            [~, FF] = fileparts(tempname);
            foo_path = fullfile(ad_local_dir, [FF, ext]);
            websave(foo_path, fileurl);
            
            % rename file when complete
            movefile(foo_path, local_path_tar)
            
            % unzip tile
            untar(local_path_tar, fullfile(ad_local_dir, fName))
            fprintf('ArcticDEM Downloaded %s: %.1fmin\n', runName, (now-time0)*24*60)
        end
        % add path to extracted DEM
        a = dir(fullfile(ad_local_dir, fName, '*dem.tif'));
        
        ad_tile(i).dem_file_name = a.name;
        ad_tile(i).dem_local_dir = fullfile(ad_local_dir, fName);
    end
else
    ad_tile(1).dem_file_name = 'WV_Greenland_Surf_210m.tif';
    ad_tile(1).dem_local_dir = [pwd,'\Data\ArcticDEM\'];
end

%% subsample Reference Ground Tracks (RGTs) for ROI
[~, laser] = arcticRGT(roi,params(1).sample_rate);

if single_track
    % select only the longest track in the ROI
    np = zeros(size(laser(1).track));
    for i = 1:length(np)
        np(i) = length(laser(1).track(i).x_RGT);
    end
    idx = find(np == max(np), 1);
    for i = 1:length(laser)
        laser(i).track = laser(i).track(idx);
    end
end

fprintf('RGTs generated for %s: %.1fmin\n', runName, (now-time0)*24*60)

%% find if point is over ice or not [TO-DO]

%% make photon cloud by sampling DEM and accounting for surface rougness,
% random probability of return, and insturment characteristics
demNoDataValue = -9999;
demFile = fullfile(ad_tile.dem_local_dir, ad_tile.dem_file_name);
% [DEM.z, R] = geotiffread(demFile);
% 
% DEM.z(DEM.z == demNoDataValue) = nan;
% DEM.x = ((R.XWorldLimits(1)+R.CellExtentInWorldX/2):R.CellExtentInWorldX:(R.XWorldLimits(2)-R.CellExtentInWorldX/2));
% DEM.y = fliplr((R.YWorldLimits(1)+R.CellExtentInWorldY/2):R.CellExtentInWorldY:(R.YWorldLimits(2)-R.CellExtentInWorldY/2))';
DEM=read_geotif_xy(demFile, roi_xlim, roi_ylim);
DEM.z=flipud(DEM.z);
DEM.y=DEM.y(end:-1:1);

c_inds = find_nearest(thinning_x,min(DEM.x)-2000):find_nearest(thinning_x,max(DEM.x)+2000);
r_inds = find_nearest(thinning_y,min(DEM.y)-2000):find_nearest(thinning_y,max(DEM.y)+2000);

% % reduce size of DEM data
% in_roi_x =  DEM.x >= (roi_xlim(1)) & DEM.x <= (roi_xlim(2));
% in_roi_y =  DEM.y >= (roi_ylim(1)) &  DEM.y <= (roi_ylim(2));
% DEM.z = DEM.z(in_roi_y,in_roi_x);
% DEM.x = DEM.x(in_roi_x);
% DEM.y = DEM.y(in_roi_y);

DEMz_orig = DEM.z;


for i = 1:length(laser)
    
    isStongBeam = any(strongBeam == i);
    for j = 1:length(laser(i).track)
        
        fracstep = laser(i).track(j).trackID/1387;
        applied_dh = (thinning(r_inds,c_inds,acquisition_step+1)-thinning(r_inds,c_inds,acquisition_step))*fracstep + thinning(r_inds,c_inds,acquisition_step);
        dhdt=interp2(thinning_x(c_inds),thinning_y(r_inds)',applied_dh,DEM.x(:)',DEM.y(:));
        %dhdt = regrid(thinning_x(c_inds),thinning_y(r_inds),applied_dh,DEM.x,DEM.y);
        
        DEM.z = DEMz_orig+dhdt;
        
        track_ind = find(tracknum == laser(i).track(j).trackID);
        mispointing_applied = mispointing(track_ind,acquisition_step);
        pe_mag_applied = pe_mag(track_ind,acquisition_step);
        pe_theta_applied = pe_theta(track_ind,acquisition_step);
               
        xy = [laser(i).track(j).X', laser(i).track(j).Y'];
        
        %%% Here we define the mispointing errors.
        cross_track_vec = make_ortholine(xy,mispointing_applied);
        pointing_error{i,j} = make_angleline(xy,pe_mag_applied,pe_theta_applied);
        
        % this section of code will run in parallel if a parpool has been opened
        if isStongBeam
            tic
            D(i).track(j) = pho_cloud_from_dem(DEM, xy+cross_track_vec+pointing_error{i,j}, params(1),optical_thickness(acquisition_step,kk)); % will sample and included dead time
            toc
            tic
            DEM1=struct('x', DEM.x,'y', DEM.y(end:-1:1),'z', DEM.z(end:-1:1,:));
            DDD=det_sim(DEM1, (xy+cross_track_vec+pointing_error{i,j})*[1; i],params(1), optical_thickness(acquisition_step,kk)+zeros(size(xy)));
            toc
            disp('done');
        else
             
            D(i).track(j) = pho_cloud_from_dem(DEM, xy+cross_track_vec+pointing_error{i,j}, params(2), optical_thickness(acquisition_step,kk)); % will sample and included dead time
        end      
    end
end


for i = 1:length(laser)  
    for j = 1:length(laser(i).track)
        
        track_ind = find(tracknum == laser(i).track(j).trackID);
        mispointing_applied = mispointing(track_ind,acquisition_step);
        pe_mag_applied = pe_mag(track_ind,acquisition_step);
        pe_theta_applied = pe_theta(track_ind,acquisition_step);
        
        D(i).track(j).x_RGT = laser(i).track(j).x_RGT(D(i).track(j).pulse_num)';
        D(i).track(j).y_RGT = laser(i).track(j).y_RGT(D(i).track(j).pulse_num)'+mispointing_applied;
        D(i).track(j).track_number = laser(i).track(j).trackID;
        %%%%%%%%%%% Here we remove knowledge of the pointing error from
        %%%%%%%%%%% the position
        D(i).track(j).x = pointing_error{i,j}(D(i).track(j).pulse_num,1);
        D(i).track(j).y = pointing_error{i,j}(D(i).track(j).pulse_num,2);
        D(i).track(j).t = (D(i).track(j).x_RGT / params(1).ground_speed);        
    end
end

clear DEM
clear laser

%% Generate ATL03 like data [TO-DO]
% Anita's ground finding code to be included here at later date
% write_ATL03_h5

%% Generate ATL06 like data

% Reorganize data for input into ATL06 ground finder
% repeat measurement index
k_rep = 1; % for future development

% making this a parfor crashes my computer
for bp = 1:3
    track = struct();
    for t = 1:length(D(1).track)
        b1 = strongBeam(bp);
        b2 = weakBeam(bp);
        
        track(t).atl03(1).BGR = ones(size(D(b1).track(t).x_RGT),'single')*params(1).NoiseRate;
        track(t).atl03(2).BGR = ones(size(D(b2).track(t).x_RGT),'single')*params(2).NoiseRate;
        
        track(t).atl03(1).beam = zeros(size(D(b1).track(t).x_RGT),'single') + b1;
        track(t).atl03(2).beam = zeros(size(D(b2).track(t).x_RGT),'single') + b2;
        
        track(t).atl03(1).track = ones(size(D(b1).track(t).x_RGT),'single')*D(b1).track(t).track_number;
        track(t).atl03(2).track = ones(size(D(b2).track(t).x_RGT),'single')*D(b2).track(t).track_number;
        
        [track(t).atl03(1).lat, track(t).atl03(1).lon]=polarstereo_inv(D(b1).track(t).x_grd,D(b1).track(t).y_grd,6378137.0,0.08181919,70,-45);
        [track(t).atl03(2).lat, track(t).atl03(2).lon]=polarstereo_inv(D(b2).track(t).x_grd,D(b2).track(t).y_grd,6378137.0,0.08181919,70,-45);
        
        track(t).atl03(1).pulse_num = single(D(b1).track(t).pulse_num);
        track(t).atl03(2).pulse_num = single(D(b2).track(t).pulse_num);
        
        track(t).atl03(1).time = D(b1).track(t).t/24/3600+91*(acquisition_step-1)+91*(track(t).atl03(2).track(1)/1387); % time in days
        track(t).atl03(2).time = D(b2).track(t).t/24/3600+91*(acquisition_step-1)+91*(track(t).atl03(2).track(1)/1387); % time in days
        
        track(t).atl03(1).x_RGT = D(b1).track(t).x_RGT;
        track(t).atl03(2).x_RGT = D(b2).track(t).x_RGT;
        
        track(t).atl03(1).seg_num=ceil(track(t).atl03(1).x_RGT/20);
        track(t).atl03(2).seg_num=ceil(track(t).atl03(2).x_RGT/20);
        
        track(t).atl03(1).y_RGT = single(D(b1).track(t).y_RGT);
        track(t).atl03(2).y_RGT = single(D(b2).track(t).y_RGT);
        
        track(t).atl03(1).h = single(D(b1).track(t).h);
        track(t).atl03(2).h = single(D(b2).track(t).h);
        
        track(t).atl03(1).h0 = single(D(b1).track(t).z);
        track(t).atl03(2).h0 = single(D(b2).track(t).z);
        
        track(t).atl03(1).detected = single(D(b1).track(t).recorded);
        track(t).atl03(2).detected = single(D(b2).track(t).recorded);
        
        track(t).atl03(1).t_ph = single(D(b1).track(t).t_ph);
        track(t).atl03(2).t_ph = single(D(b2).track(t).t_ph);
        
        % 0-noise. 1- added to allow for buffer but algorithm thinks is noise,-2-low, 3-med, 4-high
        % --- for now let ATLAS_L3a_proc_ATBD do the atl03 gorundfinding ---
        track(t).atl03(1).ph_class = zeros(size(D(b1).track(t).x_RGT), 'single');
        track(t).atl03(2).ph_class = zeros(size(D(b1).track(t).x_RGT), 'single');
    end
    
    %         x_start = track(t).atl03(2).x_RGT(1);
    %         hold off
    %         plot((track(t).atl03(1).x_RGT - x_start)/1E3, track(t).atl03(1).h, '.');
    %         hold all
    %         plot((track(t).atl03(2).x_RGT - x_start)/1E3, track(t).atl03(2).h, '.');
    %
    %         keyboard
    %
    
    % save atl03_sim beam pair data
    atl03_sim_fName = fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', ['atl03_sim2_',sprintf('%0.2d',acquisition_step)], [runName '_BP' num2str(bp) '.mat']);
    save_atl06_sim(atl03_sim_fName, track)
end

clear D
fprintf('atl03 like data generated for %s: %.1fmin\n', runName, (now-time0)*24*60)

%% apply ATL06 ground finder
% this take about 20 min
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields)
    SNR_F_table.(fields{kf})= single(h5read(snr_f_file, ['/',fields{kf}]));
end


if length(track(t)) > 0
    % if this takes up too much memeory make this loop a "for" and next a
    % "parfor"
    
    
    for bp = 1:3
        track = [];
        % simulated atl03 can be very large, load in subsets
        atl03_sim_fName = fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', ['atl03_sim2_',sprintf('%0.2d',acquisition_step)], [runName '_BP' num2str(bp) '.mat']);
        A3 = load(atl03_sim_fName, 'track');
        
        for t = 1:length(A3.track)
            
            if true
                track(t).atl06 = ATLAS_L3a_proc_ATBD(A3.track(t).atl03, params, [], SNR_F_table);
            else
                atl06_sim_fName = fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', 'atl06_sim', [runName '_BP' num2str(bp) '.mat']);
                H = load(atl06_sim_fName, 'track');
                track = H.track;
            end
            
            if length(track(t).atl06) > 0
                % add truth topo values
                track(t).atl06.h_mean_truth = nan(size(track(t).atl06.h_mean));
                track(t).atl06.h_med_truth = nan(size(track(t).atl06.h_mean));
                track(t).atl06.h_med_truth = nan(size(track(t).atl06.h_mean));
                track(t).atl06.bsnow_od = optical_thickness(acquisition_step,kk);
                for b = 1:2
                    
                    for i = 1:size(track(t).atl06.segment_id,1);
                        seg = track(t).atl06(1).segment_id(i,1);
                        idx = A3.track(t).atl03(b).seg_num == seg;
                        track(t).atl06.h_mean_truth(i,b) = mean(A3.track(t).atl03(b).h0(idx), 'omitnan');
                        track(t).atl06.h_med_truth(i,b) = median(A3.track(t).atl03(b).h0(idx), 'omitnan');
                    end
                end
            end
        end
        
        % save atl06_sim beam pair data
        atl06_sim_fName = fullfile('/','Volumes','ice2','nick','34_IceSat2','SyntheticCode_Alex','Data', ['atl06_sim_',sprintf('%0.2d',acquisition_step)], [runName '_BP' num2str(bp) '.mat']);
        save_atl06_sim(atl06_sim_fName, track)
        fprintf('atl06 beam pair %.0f complete for %s: %.1fmin\n', bp, runName, (now-time0)*24*60)
    end
    fprintf('atl06 like data generated for %s: %.1fmin\n', runName, (now-time0)*24*60)
    

end

total_time = toc;
save('Benchmark.mat','total_time');

end