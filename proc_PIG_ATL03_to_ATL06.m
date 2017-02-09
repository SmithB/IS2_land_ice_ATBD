% name these conditions
run_type='dh_clouds';


% setup the path defaults
IS_paths=IS_LI_paths('dh_clouds');

% Define the cloud fraction PDF
% optical thicknesses:
Pcloud.tau={'0.0' '1.0' '2.0' '3.0' 'Nothing'};
% probability of clouds at each thickness:
Pcloud.P=[24.6 37.5/2 37.5/2 37.9/2 37.9/2]/100;

% load the ground tracks:
load data_files/PIG_groundtracks

% choose how many repeats to simulate:
Nreps=12;
% choose how many times to run the simluation
N_runs=10;
  

if false
    [s, out]=unix(sprintf('ls %s/tau=*/rep_*/*sig*.h5',IS_paths)); out=strsplit(deblank(out));
    bad_sig_files={};
    good_sig_flag=true(length(out),1);
    for k=1:length(out)
        try
            I=h5info(out{k});
        catch
            bad_sig_files{k}=out{k};
            good_sig_flag(k)=false;
        end
    end
    bad_sig_files=bad_sig_files(~good_sig_flag);
    if ~isempty(bad_sig_files)
        for k=1:length(bad_sig_files)
            fprintf(1,'%s\n', bad_sig_files{k});
            delete(bad_sig_files{k});
        end
    end
end

load WF_est
tau_list={'0.0','1.0','2.0','3.0'};tau_list=tau_list(end:-1:1);
badlist={};
for k_tau=length(tau_list):-1:1
    top_dir=sprintf(['%s/tau=%s/'],IS_paths.ATL03,  tau_list{k_tau});
     
    d_rep=dir([top_dir,'/rep_*']);
    for k_rep=1:length(d_rep);
        rep_dir=[top_dir,'/', d_rep(k_rep).name];
        %        d_TP=[dir([rep_dir,'/*.mat']);dir([rep_dir,'/*D2.h5'])];
        d_TP=[dir([rep_dir,'/*.mat']);dir([rep_dir,'/*Track_*D2.h5'])];
        for k_TP=1:length(d_TP);
            out_dir=strrep(rep_dir,'ATL03','ATL06');
            in_file=[rep_dir,'/', d_TP(k_TP).name];
            out_file=strrep(in_file,'ATL03','ATL06');
            out_file=strrep(out_file, 'D2.h5','D3.h5');
            
            % run only on files in a short list appropriate for the current
            % subset
            if exist('top_level_params','var') && isfield(top_level_params,'files');
                [~, thebase]=fileparts(out_file);
                if ~ismember(thebase, top_level_params.files);
                    continue; 
                end
            end
            
            if ~exist(out_dir,'dir');
                mkdir(out_dir);
            end
            if exist(out_file,'file'); continue; end
            status=lockfile_tool('lock', out_file);
            if status~=0
                continue
            end
            
            fprintf(1, 'processing %s\n\tto \t%s\n', in_file, out_file);
            [D2a, PairData, params, TrackData]=read_ATLAS_h5_D2a(in_file);
            
            if exist('top_level_params','var') && isfield(top_level_params,'XR');
                for kB=1:2;
                    [xx,yy]=ll2ps(D2a(kB).lat, D2a(kB).lon);
                    D2a(kB)=index_struct(D2a(kB), xx > top_level_params.XR(1) & xx < top_level_params.XR(2) & yy > top_level_params.YR(1) & yy < top_level_params.YR(2));
                end         
            end
                
                
            if isempty(D2a);
                continue;
            end    
            
             % patch to fix the field names
            if ~isfield(params, 'N_channels');
                [params(:).N_channels]=deal(params(:).N_det);
                params=rmfield(params,'N_det');
            end
             
            for k=1:2;
                params(k).cycle=k_rep;
                params(k).RGT=median(PairData.track);
                params(k).PT=median(PairData.pair);
                params(k).GT=2*(params(k).PT-1)+k;
                params(k).orbit_number=(params(k).cycle-1)*1387+params(k).RGT;
                params(k).ATL03_sig_find=true;
            end
            [D3, dh_hist]=ATLAS_L3a_proc_ATBD(D2a,  params);
            if isempty(D3); continue; end
            
            if ~isfield(PairData,'xy');
                PairData.xy=PairData.x+1i*PairData.y;
            end
             
            write_ATL06_h5(D3, dh_hist, out_file);
            %save(out_file,'D3');
            lockfile_tool('unlock', out_file);
        end
    end
end