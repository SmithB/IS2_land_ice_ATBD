function [out_file, filename]=proc_one_GL_sim_file(filename, skip_proc)

[out_file, filename, this_cycle, this_RGT]=make_paths(filename);
if exist('skip_proc','var') && skip_proc; return; end

load WF_est
params_both=struct('RGT', this_RGT, 'WF', WF, 'cycle', this_cycle, 'pulses_per_seg', 57, 'sigma_x', 7.5, 't_dead', 3.2e-9);
params=repmat(params_both, 1, 2);
                                                                                                                     
% read in the SNR F table:                                                                                             
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};                                                       
for kf=1:length(fields)                                                                                                
    SNR_F_table.(fields{kf})=h5read('SNR_F_table_June26_2018.h5', ['/',fields{kf}]);                                   
end       

params(1).N_channels=16;
params(1).N_det=16;
params(2).N_channels=4;
params(2).N_det=4;

beam_names={'gt1l','gt1r','gt2l','gt2r','gt3l','gt3r'};

for pair=1:3
    D2a=read_ATLAS_h5_D2a(filename, pair);
    if isempty(D2a) || isempty(D2a(1).h)  || isempty(D2a(2).h)
        continue
    end
    for k=1:2
        params(k).PT=pair;
    end
    D2a=launder_ATL03(D2a);  
    
    [temp1, temp2]=ATLAS_L3a_proc_ATBD('ATL03_to_ATL06', D2a, params, [], SNR_F_table);
    
    if ~isempty(temp1)
        [D3(pair), dh_hist{pair}]=deal(temp1, temp2);
    end
end

if exist('D3','var')
    save(out_file, 'D3','dh_hist')
end


if false
    files=glob('/Volumes/ice2/ben/GreenlandSimGF/v6/*/*sigparm*.h5');
    fid=fopen('ATL06_queue.txt','w');
    for k=1:length(files)
        disp(files{k})
        [out_file, filename]=proc_one_GL_sim_file(files{k}, true);
        if ~exist(out_file,'file')
            fprintf(fid,'proc_one_GL_sim_file(''%s'');\n', files{k});
        end
    end
end



function [out_file, filename, this_cycle, this_RGT]=make_paths(filename)


filename=strrep(filename,'//','/');
filename=strrep(filename, 'V052015.0_sigparms.h5','');

temp=regexp(filename,'/TrackData_v.*_(\d+).*Track_(\d)+.*.h5','tokens');
this_cycle=str2num(temp{1}{1});
this_RGT=str2num(temp{1}{2});
[thedir, thefile, ~ ]=fileparts(filename);
[base_dir, cycle_dir]=fileparts(thedir);
out_dir=[base_dir,'/ATL06_v6/',cycle_dir];
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
out_file=[out_dir,'/', thefile,'.mat'];
