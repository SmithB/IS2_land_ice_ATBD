 function D3a=test_one_SNR_significance_val(tau, BGR, seg_center_vals, N_pulses, N_chan, roughness, N_per_pulse_0, Htot, WF)

%fprintf(1,'starting test_SNR_significance with BGR=%3.3d MHz, roughness=%3.2d\n', BGR/1e6, roughness); tic
tic
DOPLOT=false;

LS_fit_options=struct( 'Nsigma', 3, 'Hwin_min', 3, 'restrict_fit_to_initial_els', false,'SAVELOG', true);


[D2, params]=make_ATL03_data(N_pulses, N_chan, roughness, N_per_pulse_0*exp(-tau), WF, BGR, Htot, 0);
D2.seg_num=floor((D2.pulse_num/(57/2)));
params.sigma_pulse=0.68e-9;
D2.x_RGT=D2.x0;


% make a dummy structure;
temp=make_ATL03_data(57, N_chan, 0, 3, WF, BGR, Htot, 0);
selected_PE=abs(temp.h)<5;
temp.seg_num=floor((temp.pulse_num/(57/2))); temp.x_RGT=temp.x0;
[temp]=ATLAS_L3a_proc_ATBD('ATLAS_LS_fit', temp, 57/2,selected_PE, 100, params, ...
    struct('seg_count', 1,'N_final', 0,'w_surface_window_initial', 10,'N_seg_pulses', 58), LS_fit_options);
temp.N_final=NaN;
ff=fieldnames(temp);
for kf=1:length(ff); temp.(ff{kf})=NaN; end
D3a=repmat(temp, [length(seg_center_vals), 1]);


for k0=1:length(seg_center_vals)
    this_seg=seg_center_vals(k0);
    x_seg_ctr=(this_seg*57*.7)/2;
    % hack to speed up the computation:
    temp=range(find(D2.seg_num==this_seg-1 | D2.seg_num==this_seg));
    D2sub=index_struct(D2, temp(1):temp(2) );
    D2sub.ph_class=zeros(size(D2sub.h));
    [selected_PE, signal_selection_source, signal_selection_status]=ATLAS_L3a_proc_ATBD('backup_signal_finding_strategy', D2sub, D2, this_seg, 10);
    
    if sum(selected_PE) > 10
        [D3a(k0), r, els, LOG]=ATLAS_L3a_proc_ATBD('ATLAS_LS_fit', D2sub, x_seg_ctr,selected_PE, diff(range(D2sub.h(selected_PE))), params, ...
            struct('seg_count', k0,'N_final', 0,'w_surface_window_initial', Htot, 'N_seg_pulses', 57), LS_fit_options);
        D3a(k0).N_final=sum(abs(r)<D3a(k0).w_surface_window_final<2);
        %if abs(D3a(k0).h_mean) > 2 ; DOPLOT=true; else; DOPLOT=false; end 
        if DOPLOT
            figure(k0); clf; plot(D2sub.x_RGT, D2sub.h,'r.');
            hold on; plot(D2sub.x_RGT(selected_PE), D2sub.h(selected_PE),'bo');
            plot_segs(x_seg_ctr, D3a(k0).h_mean, D3a(k0).dh_fit_dx, 40, 'r');
            D3a(k0)
            disp('');
        end
    end
end

D3a=flatten_struct(D3a);
fprintf(1,'finished test_SNR_significance with BGR=%3.1f MHz, tau=%2.1f, roughness=%2.1f in %2.1f seconds\n', BGR/1e6, tau, roughness, toc);
