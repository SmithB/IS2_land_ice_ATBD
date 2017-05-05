function  [ATL06_data, BGR_vals, Hwin_vals]=test_ATL06_noise_convergence(kBG)

% what is the false-positive rate for a given SNR and window size?

% Generate data spanning +- 100 m vertical 
% Subset to smaller window (for a range of Hwin values)
% run the guts of ATL06
% Calculate the SNR


N_MC=1e5;
load WF_est;   

BGR_vals=logspace(5.5, 7, 20);
Hwin_vals=[2 3 4 8 10 12 16 20 30 40 80 120 160 200]; 
 
if false
    fid=fopen('ATL06_BG_queue.txt','w');
    for k=1:length(BGR_vals)
        fprintf(fid,'test_ATL06_noise_convergence(%d)\n', k);
    end
    fclose(fid);
end

% parfor kBG=1:length(BGR_vals)
%     ATL06_data{kBG}=test_one_convergence_model(Hwin_vals, BGR_vals(kBG), N_MC, WF);      
% end
for k=1:ceil(N_MC/1000)
    ATL06_data_temp{k}=test_one_convergence_model(Hwin_vals, BGR_vals(kBG), WF);
end
ATL06_data=flatten_struct(cat(1, ATL06_data_temp{:}));

save(sprintf('ATL06_convergence_output/%d.mat',kBG),'ATL06_data','BGR_vals','Hwin_vals','kBG');


% build the look-up table;
if false
    %SNR_F_table is N_SNR x N_Hwin x N_BGR
    
    SNR_F_table.SNR=linspace(-10, 20, 301);
    SNR_F_table.BGR=logspace(5.5, 7, 20);
    SNR_F_table.H_win_initial=[2 3 4 8 10 12 16 20 30 40 80 120 160 200]; 
 
    for kB=1:20
        thefile=sprintf('ATL06_convergence_output/%d.mat',kB);
        if ~exist(thefile,'file')
            continue
        end
        load(thefile);
        for kH=1:size(ATL06_data.h_mean,2)
            for kS=1:length(SNR_F_table.SNR)
                SNR_F_table.F(kS, kH, kB)=mean(isfinite(ATL06_data.h_mean(:, kH)) & ATL06_data.SNR(:, kH) > SNR_F_table.SNR(kS)); 
            end
        end
    end 
    
    if false  % DON't do this until you are ready to overwrite the existing table
        h5create('SNR_F_table.h5','/P_NoiseOnly', [301, 14, 20],'datatype','double')
        h5write('SNR_F_table.h5','/P_NoiseOnly', SNR_F_table.F)
        h5create('SNR_F_table.h5','/SNR', 301 ,'datatype','double')
        h5write('SNR_F_table.h5','/SNR', SNR_F_table.SNR)
        
        h5create('SNR_F_table.h5','/BGR', 20,'datatype','double')
        h5write('SNR_F_table.h5','/BGR', SNR_F_table.BGR)
        h5create('SNR_F_table.h5','/W_surface_window_initial', 14,'datatype','double')
        h5write('SNR_F_table.h5','/W_surface_window_initial', SNR_F_table.Hwin_initial)
    end
end


%------------------------------------------------------------------------
function D3=test_one_convergence_model(Hwin_vals, BGR, WF)
  tic
N_pulses=57*1000;
N_chan=4;
Htot=200;
seg_center_vals=3:(N_pulses/57-3);

DOPLOT=false;

LS_fit_options=struct( 'Nsigma', 3, 'Hwin_min', 3);
 

[D2, params]=make_ATL03_data(N_pulses, N_chan, 0, 0, WF, BGR, Htot, 0);
D2.seg_num=floor((D2.pulse_num/(57/2)));
params.sigma_pulse=0.68e-9;
D2.x_RGT=D2.x0;
clear D3a;

% make a dummy structure;
temp=make_ATL03_data(57, N_chan, 0, 3, WF, BGR, Htot, 0);
selected_PE=abs(temp.h)<5;
temp.seg_num=floor((temp.pulse_num/(57/2))); temp.x_RGT=temp.x0;
[temp]=ATLAS_L3a_proc_ATBD('ATLAS_LS_fit', temp, 57/2,selected_PE, 100, params, ...
    struct('seg_count', 1,'N_final', 0,'w_surface_window_initial', 10), LS_fit_options);
temp.N_final=NaN;
ff=fieldnames(temp);
for kf=1:length(ff); temp.(ff{kf})=NaN; end
D3a=repmat(temp, length(seg_center_vals), length(Hwin_vals));

for k0=1:length(seg_center_vals)
    this_seg=seg_center_vals(k0);
    x_seg_ctr=(this_seg*57*.7)/2;
    D2sub=index_struct(D2, D2.seg_num==this_seg-1 | D2.seg_num==this_seg);
    D2sub.ph_class=zeros(size(D2sub.h));
    
    for kH=1:length(Hwin_vals)
        if Hwin_vals(kH) > 10
            % first, select PE using the backup signal finding strategy:
            [selected_PE, ~, ~]=ATLAS_L3a_proc_ATBD('backup_signal_finding_strategy', D2sub, D2, this_seg, 10);
%             if sum(selected_PE) > 10;
%                 disp('hup!');
%             end
        else
            selected_PE=abs(D2sub.h) < Hwin_vals(kH)/2;
        end
        if sum(selected_PE) > 10
            [D3a(k0, kH), r, els]=ATLAS_L3a_proc_ATBD('ATLAS_LS_fit', D2sub, x_seg_ctr,selected_PE, 100, params, ...
                struct('seg_count', k0,'N_final', 0,'w_surface_window_initial', diff(range(D2sub.h(selected_PE)))), LS_fit_options);
            D3a(k0, kH).N_final=sum(abs(r)<D3a(k0, kH).w_surface_window_final<2);
        end
    end
end

f=fieldnames(D3a);
for kf=1:length(f);
    this=f{kf};
    for kH=1:length(Hwin_vals);
        D3.(this)(:, kH)=[D3a(:, kH).(this)];
    end
end
toc
    



