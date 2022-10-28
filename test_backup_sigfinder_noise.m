function  [BSD_temp]=test_backup_sigfinder_noise(kBG, N_MC)

if ~exist('kBG','var') && ~exist('N_MC','var'); 
    make_queue
    return
end

% what is the false-positive rate for a given SNR and window size?

% Generate data spanning +- 100 m vertical
% Subset to smaller window (for a range of Hwin values)
% run the guts of ATL06
% Calculate the SNR

load WF_est;

BGR_vals=logspace(log10(7e5), 7, 20);
Hwin_vals=[3 4 5 7.5 10 15 20 30 40 80 120 160 200];

% parfor kBG=1:length(BGR_vals)
%     ATL06_data{kBG}=test_one_convergence_model(Hwin_vals, BGR_vals(kBG), N_MC, WF);
% end
for k=1:ceil(N_MC/1000)
    BSD_temp(k)=test_one_convergence_model(Hwin_vals, BGR_vals(kBG), WF);
end

BSD_temp=flatten_struct(BSD_temp);

out_file=sprintf('backup_sigfinder_test/test_file%d.mat', kBG);
save(out_file, 'BGR_vals', 'Hwin_vals','BSD_temp','N_MC');


%------------------------------------------------------------------------
function D_out=test_one_convergence_model(Hwin_vals, BGR, WF)
tic
N_pulses=57*1000;
N_chan=4;
Htot=200;
seg_center_vals=3:2:(N_pulses/(57/2)-3);

[D2, params]=make_ATL03_data(N_pulses, N_chan, 0, 0, WF, BGR, Htot, 0);
D2.seg_num=floor((D2.pulse_num/(57/2)));
params.sigma_pulse=0.68e-9;
D2.x_RGT=D2.x0;
clear D3a;

D2=index_struct(D2, true(size(D2.h)),{'h','seg_num','x_RGT'});

[D_out.N_HW, D_out.N_tot, D_out.HW]=deal(NaN(length(seg_center_vals), length(Hwin_vals)));

for k0=1:length(seg_center_vals)
    this_seg=seg_center_vals(k0);
    D2sub=index_struct(D2, D2.seg_num==this_seg-1 | D2.seg_num==this_seg);
    D2sub.ph_class=zeros(size(D2sub.h));
    
    D2all_seg=index_struct(D2,ismember(D2.seg_num, this_seg+[-2:2]));
    
    for kH=1:length(Hwin_vals)
        D2sub2=index_struct(D2sub, abs(D2sub.h) < Hwin_vals(kH)/2);
        D2_all_sub=index_struct(D2all_seg, abs(D2all_seg.h)<Hwin_vals(kH)/2);
        %  select PE using the backup signal finding strategy:
        [selected_PE, ~, ~]=ATLAS_L3a_proc_ATBD('backup_signal_finding_strategy', D2sub2, D2_all_sub, this_seg, 10);
        D_out.N_HW(k0, kH)=sum(selected_PE);
        D_out.N_tot(k0, kH)=sum(~selected_PE);
        if D_out.N_HW(k0, kH)>0
            D_out.HW(k0, kH)=diff(range(D2sub2.h(selected_PE)));
        end
    end
end

%toc
%--------------------------------
function make_queue

% make the queue:

BGR_vals=logspace(log10(7e5), 7, 20);
Hwin_vals=[3 4 5 7.5 10 15 20 30 40 80 120 160 200];



fid=fopen('backup_sigfinder_queue.txt','w');
for k=1:length(BGR_vals)
    % want at least 2500 examples that have 10 photons in the smallest
    % Hwin
    this_N_MC=max(1e5, 2500./(1-poisson_p_table(BGR_vals(k)/1.5e8*58*10, 10)));
    if this_N_MC <= 1e5
        fprintf(fid,'test_backup_sigfinder_noise(%d, %d); test_AT_fit_noise(%d, 1e5)\n', k, round(this_N_MC), k);
    end
end
fclose(fid);



