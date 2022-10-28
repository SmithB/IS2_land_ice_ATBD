
N_reps=12;

PairTrackList=load('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_pairtrack_list');
 

if false

% make_PIG_fakedata_queue;
N_reps=12;
rep_list=1:N_reps;
tau_vals=[0  1  2  3  4];
TP_list=1:32;
fid=fopen('PIG_ATL06_queue.txt','w');

for k_rep=1:length(rep_list)
    for k_TP=1:length(TP_list)       
        for k_tau=1:length(tau_vals)
            fprintf(fid,'tau_vals=%f; rep_list=%d; TP_list=%d; Generate_PIG_fake_data;\n', tau_vals(k_tau), rep_list(k_rep), TP_list(k_TP));
        end
    end
end
fclose(fid);

end
