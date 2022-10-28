load WF_est;
%load prelaunch_wf_est_tep
sigma_0=0.5*diff(wf_percentile(WF.t, WF.p, [0.16 0.84]));
W_pulse=sigma_0*1.5e8+[0 0.025 0.05 0.1 0.2 0.4 0.8 1.6 ];

HW_vals=3+([0 0.1 0.5 1, 2, 6, 10]);

SNR_vals=[1 2 10 100];

% the minimum SNR is 10 PE /( 10 MHz * 20 m /1.5e8 m/s * 57 pulses) =
% 1.3

% need a quick test of robust_peak_width_from_hist

%SNR_vals=500;
for kWp=1:length(W_pulse)
    for kHW=1:length(HW_vals)
        for kSNR=1:length(SNR_vals)
            [dM(kWp,  kHW, kSNR), dCtr(kWp, kHW, kSNR), syn_wf(kWp, kHW, kSNR)]=correct_for_TX_shape(WF.t, WF.p, HW_vals(kHW)/1.5e8,  W_pulse(kWp)/1.5e8, SNR_vals(kSNR));            
        end
    end
end
TX_shape_corr_table.W_pulse=W_pulse;
TX_shape_corr_table.SNR=SNR_vals;
TX_shape_corr_table.HW=HW_vals;
TX_shape_corr_table.delta_med=dM;
TX_shape_corr_table.delta_mean=dCtr;

Tstr=strsplit(datestr(now));

 
 %h5_file=sprintf('/Volumes/ice2/ben/sdt/test_data/TX_shape_corr_test_data_TEP_%s.h5', Tstr{1})
 h5_file=sprintf('/Volumes/ice2/ben/sdt/test_data/TX_shape_corr_test_data_ATBD_WF_%s.h5', Tstr{1})
 
 f=fieldnames(TX_shape_corr_table);
 for kf=1:length(f)
     h5create(h5_file,['/', f{kf}], size(TX_shape_corr_table.(f{kf})),'datatype','double');
     h5write(h5_file,['/', f{kf}], TX_shape_corr_table.(f{kf}));
 end
h5create(h5_file,'/WF/t', size(WF.t), 'datatype','double');
h5write(h5_file,'/WF/t', WF.t);
h5create(h5_file,'/WF/p', size(WF.p), 'datatype','double');
h5write(h5_file,'/WF/p', WF.p);

syn_wf_fields={'t','P', 'mask'};
for kWp=1:length(W_pulse)
    for kHW=1:length(HW_vals)
        for kSNR=1:length(SNR_vals)
            group_name=sprintf('/syn_wf/iWp=%d_iHW=%d_iSNR=%d', kWp, kHW, kSNR);
            for kf=1:length(syn_wf_fields)
                ds_name=[group_name,'/',syn_wf_fields{kf}];
                h5create(h5_file, ds_name, size(syn_wf(kWp, kHW, kSNR).(syn_wf_fields{kf})),'datatype','double'); 
                h5write(h5_file, ds_name,  double(syn_wf(kWp, kHW, kSNR).(syn_wf_fields{kf})));
            end
        end
    end
end

            




