function [D2a, PairData, params, TrackData]=read_ATLAS_h5_D2a(h5_file, skip_sig)

if ~exist('skip_sig','var')
    skip_sig=false;
end

this_dir=fileparts(h5_file);
sig_file=dir([h5_file,'V*.h5']);
[~, ind]=sort([sig_file.datenum]);
if isempty(ind) && ~skip_sig
    disp(['no sig file for ' , h5_file]); 
    D2a=[]; PairData=[]; params=[]; TrackData=[]; return;
end
if ~skip_sig
    sig_file=sig_file(ind(end));
end
try
    I=h5info([this_dir,'/',sig_file.name]);
catch
    if ~skip_sig
        disp(['bad sig file for ' , h5_file]);
        D2a=[]; PairData=[]; params=[]; return;
    end
end

top_groups={'/D2_weak','/D2_strong'};
for kG=1:length(top_groups);
    I=h5info(h5_file,top_groups{kG});
    for kD=1:length(I.Datasets);
        D2a(kG).(I.Datasets(kD).Name)=h5read(h5_file,[top_groups{kG},'/', I.Datasets(kD).Name]);
    end
end
  
sig_groups={ '/channelD2_weak/photon','/channelD2_strong/photon'};
if ~skip_sig
    for kG=1:length(sig_groups);
        temp=h5read([this_dir,'/',sig_file.name],[sig_groups{kG},'/ph_class']);
        D2a(kG).ph_class=zeros(size(D2a(kG).h));
        D2a(kG).ph_class(isfinite(D2a(kG).h))=temp;
    end
end

I=h5info(h5_file,'/PairData');

for kD=1:length(I.Datasets);
   PairData.(I.Datasets(kD).Name)=h5read(h5_file,['/PairData/', I.Datasets(kD).Name]);
end

WF.p=h5read(h5_file,'/WF/p');
WF.t=h5read(h5_file,'/WF/t');

params_L=struct('N_per_pulse', 12, 't_dead', 3.2e-9, 'sigma_x', 7.2,'sigma_pulse', 1.6e-9,'c', 3e8, 'N_det', 16, 'NoiseRate', 1e7,'H_window', 80, 'WF', WF,'DEBUG', false);
params_R=params_L; params_R.N_per_pulse=3; params_R.N_det=4;
params=[params_R params_L];

try
    I=h5info(h5_file,'/TrackData');
catch
    I=[];
end
if ~isempty(I); 
    for kD=1:length(I.Datasets);
        TrackData.(I.Datasets(kD).Name)=h5read(h5_file,['/TrackData/', I.Datasets(kD).Name]);
    end
end







