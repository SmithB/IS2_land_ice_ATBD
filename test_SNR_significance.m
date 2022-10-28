% what is the false-positive rate for a given SNR_significance threshold?
% Need to test as a fn of signal strength, BGR, sigma

% Generate +- 50 m of data, 
% Center a ground window on +- 10 m
% run the guts of ATL06
% evaluate the results

 
%tau_values=4:-.25:0;
%BGR_values=[0.5:0.5:10]*1e6;
roughness_vals=[0 1 2];
tau_values=[4:-.25:0];
BGR_values=[0.25:0.25:10]*1e6;

N_pulses=57*1606;
N_chan=4;
N_per_pulse_0=3;
Htot=200;

 
load WF_est;

seg_center_vals=3:(N_pulses/57-3);

clear D3 D3a D3b;
for kR=1:length(roughness_vals)
    for kTau=1:length(tau_values)
        parfor kB=1:length(BGR_values) 
            tau=tau_values(kTau);
            BGR=BGR_values(kB);
            roughness=roughness_vals(kR);
            D3b{kTau, kB, kR}=test_one_SNR_significance_val(tau, BGR, seg_center_vals, N_pulses, N_chan, roughness, N_per_pulse_0, Htot, WF);
        end
    end
end



clear Fgood F95 F99 Fbad
for kR=1:numel(roughness_vals)
    for kTau=1:numel(tau_values)
        for kB=1:numel(BGR_values)
            D3b{kTau, kB, kR}.SNR_significance=...
                interpn(SNR_F_table.SNR, SNR_F_table.W_surface_window_initial, SNR_F_table.BGR, SNR_F_table.P_NoiseOnly, ...
                max(-10, min(20, D3b{kTau, kB, kR}.SNR)), ...
                max(3, min(max(SNR_F_table.W_surface_window_initial), Htot))*ones(size(D3b{kTau, kB, kR}.SNR)), ...
                max(1e5, min(1e7, D3b{kTau, kB, kR}.BGR)));
            good=abs(D3b{kTau, kB, kR}.h_mean)<1 & isfinite(D3b{kTau, kB, kR}.h_mean)& isfinite(D3b{kTau, kB, kR}.h_mean)~=0  & abs(D3b{kTau, kB, kR}.dh_fit_dx) < 1/20  ;
            Fgood(kTau, kB, kR)=mean(good);
            F95(kTau, kB, kR)=mean(D3b{kTau, kB, kR}.SNR_significance < 0.05 & good);
            F99(kTau, kB, kR)=mean(D3b{kTau, kB, kR}.SNR_significance < 0.01 & good);
            FgoodRejected95(kTau, kB, kR)=mean(D3b{kTau, kB, kR}.SNR_significance >0.05 & good);
            FgoodRejected99(kTau, kB, kR)=mean(D3b{kTau, kB, kR}.SNR_significance >0.01 & good);
            bad=abs(D3b{kTau, kB, kR}.h_mean)>1 & isfinite(D3b{kTau, kB, kR}.h_mean);
            Fbad(kTau, kB, kR)=mean(bad);
            FBadRejected95(kTau, kB, kR)=mean(D3b{kTau, kB, kR}.SNR_significance>0.05 & bad);
            FBadRejected99(kTau, kB, kR)=mean(D3b{kTau, kB, kR}.SNR_significance>0.01 & bad);
            FBadAccepted95(kTau, kB, kR)=mean(D3b{kTau, kB, kR}.SNR_significance<0.05 & bad);
            FBadAccepted99(kTau, kB, kR)=mean(D3b{kTau, kB, kR}.SNR_significance<0.01 & bad);
        end
    end
end






% map of the fraction of good segments rejected at each confidence level
figure(1); clf
hax=cheek_by_jowl(3,3, [0.15 0.15 0.7 0.7]); 
for k=1:3
    axes(hax(k,1)); imagesc(BGR_values/1e6, tau_values, Fgood(:,:,k)); caxis([0.0 1]); 
    axes(hax(k,2)); imagesc(BGR_values/1e6, tau_values, FgoodRejected95(:,:,k)); caxis([0 0.3]);
    axes(hax(k,3)); imagesc(BGR_values/1e6, tau_values, FgoodRejected99(:,:,k)); caxis([0 0.3]);
end; figure(gcf); 
axes(hax(1,1)); title('F good'); hb=colorbar('north');
axes(hax(1,2)); title('F good & rejected @95%'); hb=colorbar('north');
axes(hax(1,3)); title('F good & rejected @99%'); hb=colorbar('north');
axes(hax(2,1)); ylabel('optical thickness'); 
axes(hax(3,2)); xlabel('BGR, MHz');
colormap(jet*.6+.4)

% map of the fraction of bogus segments accepted at each confidence level
figure(2); clf
hax=cheek_by_jowl(3,3, [0.15 0.15 0.7 0.7]); 
for k=1:3
    axes(hax(k,1)); imagesc(BGR_values/1e6, tau_values, Fbad(:,:,k)); caxis([0.0 1]); 
    axes(hax(k,2)); imagesc(BGR_values/1e6, tau_values, FBadAccepted95(:,:,k)); caxis([0 0.3]);
    axes(hax(k,3)); imagesc(BGR_values/1e6, tau_values, FBadAccepted99(:,:,k)); caxis([0 0.3]);
end; figure(gcf);
axes(hax(1,1)); title('F bad'); hb=colorbar('north');
axes(hax(1,2)); title('F bad & accepted @95%'); hb=colorbar('north');
axes(hax(1,3)); title('F bad & accepted @99%'); hb=colorbar('north');
axes(hax(2,1)); ylabel('optical thickness'); 
axes(hax(3,2)); xlabel('BGR, MHz');
colormap(jet*.6+.4)

figure(3);  clf
hax=cheek_by_jowl(3,3, [0.15 0.15 0.7 0.7]); 
for k=1:3
    axes(hax(k,1)); imagesc(BGR_values/1e6, tau_values, Fbad(:,:,k)); caxis([0. 1]); 
    axes(hax(k,2)); imagesc(BGR_values/1e6, tau_values, FBadRejected95(:,:,k)); caxis([0 1]);
    axes(hax(k,3)); imagesc(BGR_values/1e6, tau_values, FBadRejected99(:,:,k)); caxis([0 1]);
end; figure(gcf);

axes(hax(1,1)); title('F bad'); hb=colorbar('north');
axes(hax(1,2)); title('F bad & rejected @95%'); hb=colorbar('north');
axes(hax(1,3)); title('F bad & rejected @99%'); hb=colorbar('north');
axes(hax(2,1)); ylabel('optical thickness'); 
axes(hax(3,2)); xlabel('BGR, MHz');
colormap(jet*.6+.4)

