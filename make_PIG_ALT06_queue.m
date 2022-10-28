
 
PairTrackCombos=load('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_pairtrack_list');
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


% make_PIG_fakedata_queue;
N_reps=12;
rep_list=N_reps:-1:1;
tau_vals=[0  1  2  3  4];
TP_list=1:31;
fid=fopen('PIG_ATL06_queue.txt','w');
%  fid=1;
for k_rep=1:length(rep_list)
    for k_TP=1:length(TP_list)
        for k_tau=1:length(tau_vals)
            top_dir=sprintf('%s_8.00MHz/tau=%1.1f/',IS_paths.ATL03,  tau_vals(k_tau));
            rep_dir=sprintf('%s/rep_%d/', top_dir, rep_list(k_rep));
            
            in_file=sprintf('%s/Track_%d-Pair_%d_D2.h5', rep_dir, PairTrackCombos.track(k_TP),  PairTrackCombos.pair(k_TP));
            if ~exist(in_file,'file')
                fprintf(1,'    %s doesn''t exist\n', in_file);
                continue; 
            end
            out_dir=strrep(rep_dir,'ATL03','ATL06');
            out_file=strrep(in_file,'ATL03','ATL06');
            out_file=strrep(out_file, 'D2.h5','D3.h5');
            if exist(out_file,'file')
                fprintf(1,'   --- %s exists\n', out_file)
                continue; 
            end
            fprintf(fid,'proc_one_ATL03_to_ATL06(''%s'');\n', in_file);
            
            
            %fprintf(fid,'tau_list={''%1.1f''}; rep_list=%d; TP_list=%d; proc_PIG_ATL03_to_ATL06;\n', tau_vals(k_tau), rep_list(k_rep), TP_list(k_TP));
        end
    end
end
fclose(fid);

 