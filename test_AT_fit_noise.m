function test_AT_fit_noise(kBG, N_MC_AT)

backup_sigfinder_file=sprintf('backup_sigfinder_test/test_file%d.mat', kBG);
save_file=sprintf('AT_fit_noise_test_June26_2018/test_file%d.mat', kBG);
load(backup_sigfinder_file)


BGR=BGR_vals(kBG);

LS_fit_options=struct( 'Nsigma', 3, 'Hwin_min', 3,'restrict_fit_to_initial_els', false);
params=struct('sigma_x', 7.5,'sigma_pulse', 0.68e-9); 

x_seg_ctr=20;

[out.N, out.SNR, out.HW_final, out.P0]=deal(NaN([length(Hwin_vals), N_MC_AT]));
for k_HW=1:length(Hwin_vals)
    
    good=isfinite(BSD_temp.HW(:, k_HW)) & BSD_temp.N_HW(:, k_HW)>10;
    if ~any(good)
        continue
    end
    this_P0=mean(good);
    HW_samp=BSD_temp.HW(good, k_HW);
    N_HW_samp=BSD_temp.N_HW(good, k_HW);
    % Note: N_tot is the number of unselected photons(face-palm)
    N_tot_samp=BSD_temp.N_tot(good, k_HW);
    for k_MC=1:N_MC_AT
        this=ceil(rand(1)*length(N_HW_samp));
        N_unsel=N_tot_samp(this);
        N_sel=N_HW_samp(this);
        h_in_win=((rand(N_sel,1)-0.5))*HW_samp(this);
        
        HW_outside=Hwin_vals(k_HW)-HW_samp(this);
        if HW_outside > 0 && N_unsel>0
            h_out_win=[HW_outside/2*rand(ceil(N_unsel/2),1)+HW_samp(this)/2;
                -(HW_outside/2*rand(ceil(N_unsel/2),1)+HW_samp(this)/2)];
        else
            h_out_win=[];
        end
        selected_PE=[true(size(h_in_win)); false(size(h_out_win))];
        h=[h_in_win; h_out_win];
         
        
        D2a=struct('h', h,'x_RGT', rand(size(h))*40,'BGR', zeros(size(h))+BGR);
        [D3a, r, els]=ATLAS_L3a_proc_ATBD('ATLAS_LS_fit', D2a, x_seg_ctr,selected_PE, diff(range(D2a.h(selected_PE))), params, ...
            struct('seg_count', 1,'N_final', 0,'w_surface_window_initial', HW_samp(this),'N_seg_pulses', 58), LS_fit_options);
        N=sum(abs(r)<D3a.w_surface_window_final/2);
        N_noise=D3a.w_surface_window_final*BGR/1.5e8*57;
        
        out.N(k_HW, k_MC)=N;
        out.SNR(k_HW, k_MC)= (N-N_noise)/N_noise;
        out.HW_final(k_HW, k_MC)= D3a.w_surface_window_final;
        out.P0(k_HW, k_MC)= this_P0;
    end
end

save(save_file,'out','BGR', 'Hwin_vals')

% make the queue:
if false
    fid=fopen('queue_for_AT_fit.txt','w');
    for kk=1:20;
        fprintf(fid, 'test_AT_fit_noise(%d, 5e4);\n', kk);
        
    end
    fclose(fid);
end
