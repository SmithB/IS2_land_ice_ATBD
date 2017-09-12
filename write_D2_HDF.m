function write_D2_HDF(in, out_dir)

if isstruct(in)
     [D2a, params, TrackData, PairData, file]=deal(in.D2a, in.params, in.TrackData, in.PairData, in.file);
else
    file=in;
    load(file);
end
[~,fname, ~]=fileparts(file);

if isempty(findstr(fname,'_D2'))
    h5_file=[out_dir,'/', fname,'_D2.h5'];
else
    h5_file=[out_dir,'/', fname,'.h5'];
end
if exist(h5_file,'file');
    return; %delete(h5_file);
end
disp(h5_file);
if ~exist(out_dir,'dir');
    mkdir(out_dir); 
end


% OLD HACK- for files that had a complex x_rgt 
% write out the D2 data
%fields=fieldnames(D2a);
% run through the fieldnames, convert complex variable to x and y
% for kf=1:length(fields);
%     if fields{kf}(1)~='x' 
%         continue
%     end
%     for kb=1:2
%         temp=D2a(kb).(fields{kf});
%         D2a(kb).(fields{kf})=real(temp);         
%         D2a(kb).(['y',fields{kf}(2:end)])=imag(temp);
%     end
% end

fields=fieldnames(D2a);
for k=1:length(fields); temp_f{k}=fields{k}(2:end); end; [~, ind]=sort(temp_f); fields=fields(ind);

param_fields=fieldnames(params); 
param_fields=param_fields(~ismember(param_fields,{'WF','DEBUG','c'}));

for k=1:length(fields);
    this_field=['/D2_weak/',fields{k}]; 
    temp=double(D2a(1).(fields{k}));
    h5create(h5_file, this_field, size(temp),'ChunkSize', [min(size(temp,1), 1e4), 1], 'Datatype','double','Deflate', 9);
    h5write(h5_file, this_field,  temp);
end

for k=1:length(param_fields);
    h5writeatt(h5_file,'/D2_weak/', param_fields{k}, params(1).(param_fields{k}));
end


for k=1:length(fields);
    this_field=['/D2_strong/',fields{k}]; 
    temp=double(D2a(2).(fields{k}));
    h5create(h5_file, this_field, size(temp),'ChunkSize', [min(size(temp,1), 1e4), 1], 'Datatype','double','Deflate', 9);
    h5write(h5_file, this_field,  temp);
end


for k=1:length(param_fields);
    h5writeatt(h5_file,'/D2_strong/', param_fields{k}, params(2).(param_fields{k}));
end

h5create(h5_file, '/WF/t', size(params(1).WF.t), 'Datatype', 'double');
h5write(h5_file,'/WF/t', params(1).WF.t);
h5create(h5_file, '/WF/p', size(params(1).WF.t), 'Datatype', 'double');
h5write(h5_file,'/WF/p', params(1).WF.p);


% write out the track data
fields=fieldnames(PairData);
for kf=1:length(fields);
    if fields{kf}(1)~='x';
        continue
    end         
    temp=PairData(1).(fields{kf});
    if strcmp(fields{kf},'xy');
        fields{kf}='x';
        PairData=rmfield(PairData,'xy');
    end
    PairData(1).(fields{kf})=real(temp);
    PairData(1).(['y',fields{kf}(2:end)])=imag(temp);
end

fields=fieldnames(PairData);
for kf=1:length(fields);
    h5create(h5_file,['/PairData/', fields{kf}], size(PairData(1).t), 'Datatype','double');
    h5write(h5_file,['/PairData/', fields{kf}] , PairData(1).(fields{kf}));
end


% write out the track data
fields=fieldnames(TrackData);
for kf=1:length(fields);
    if fields{kf}(1)~='x';
        continue
    end         
    temp=TrackData(1).(fields{kf});
    if strcmp(fields{kf},'xy');
        fields{kf}='x';
        TrackData=rmfield(TrackData,'xy');
    end
    TrackData(1).(fields{kf})=real(temp);
    TrackData(1).(['y',fields{kf}(2:end)])=imag(temp);
end

fields=fieldnames(TrackData);
for kf=1:length(fields);
    h5create(h5_file,['/TrackData/', fields{kf}], size(TrackData(1).t), 'Datatype','double');
    h5write(h5_file,['/TrackData/', fields{kf}] , TrackData(1).(fields{kf}));
end



