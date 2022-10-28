function proc_one_ATL03_to_ATL06(in_file, skip_sigfinder)


% read in the SNR F table:
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields)
    SNR_F_table.(fields{kf})=h5read('SNR_F_table.h5', ['/',fields{kf}]);
end

% load the ground tracks:
load data_files/PIG_groundtracks

load WF_est
 
[rep_dir, ~]=fileparts(in_file);

out_dir=strrep(rep_dir,'ATL03','ATL06');
out_file=strrep(in_file,'ATL03','ATL06');
out_file=strrep(out_file, 'D2.h5','D3.h5');

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
fprintf(1, 'processing %s\n\tto \t%s\n', in_file, out_file);
if ~exist('skip_sigfinder','var') || ~skip_sigfinder
    sig_file=dir([in_file,'*sigparms.h5']);
    sigfinder_failed=false;
    if isempty(sig_file)
        sigfinder_failed=true;
    end
    if sigfinder_failed
        return
    end
    
    [D2a, PairData, params, TrackData]=read_ATLAS_h5_D2a(in_file);
else
    [D2a, PairData, params, TrackData]=read_ATLAS_h5_D2a(in_file, true);
end


if isempty(D2a)
    return
end

% patch to fix the field names
if ~isfield(params, 'N_channels');
    [params(:).N_channels]=deal(params(:).N_det);
    params=rmfield(params,'N_det');
end


temp=regexp(out_file, 'rep_(\d+)/','tokens');
this_rep=str2double(temp{1}{1});

for k=1:2
    params(k).cycle=this_rep;
    params(k).RGT=median(PairData.track);
    params(k).PT=median(PairData.pair);
    params(k).GT=2*(params(k).PT-1)+k;
    params(k).orbit_number=(params(k).cycle-1)*1387+params(k).RGT;
    params(k).ATL03_sig_find=true;
end
[D3, dh_hist]=ATLAS_L3a_proc_ATBD(D2a,  params, [], SNR_F_table);
if isempty(D3); return; end

if ~isfield(PairData,'xy');
    PairData.xy=PairData.x+1i*PairData.y;
end

write_ATL06_h5(D3, dh_hist, out_file);



