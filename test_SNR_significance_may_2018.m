% what is the false-positive rate for a given SNR_significance threshold?
% Need to test as a fn of signal strength, BGR, sigma

% Generate +- 50 m of data, 
% Center a ground window on +- 10 m
% run the guts of ATL06
% evaluate the results

this_cmap=my_rgb_cpt(128).*repmat(linspace(0.5, 1, 128)', [1, 3]);
                                                                                                                 
% read in the SNR F table:                                                                                             
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};                                                       
for kf=1:length(fields)                                                                                                
    SNR_F_table.(fields{kf})=h5read('SNR_F_table_June26_2018.h5', ['/',fields{kf}]);                                   
end 


REGEN_DATA=false; 
%tau_values=4:-.25:0;
%BGR_values=[0.5:0.5:10]*1e6;
roughness_vals=[0 1 2];
tau_values=[4:-.25:0];
BGR_values=[0.25:0.25:10]*1e6;

N_pulses=57*1606;
N_chan=4;
N_per_pulse_0=3;
Htot=200;

if REGEN_DATA
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
    save test_SNR_significance_output_may_2018.mat
else
    if ~exist('roughnes_vals','var')
        load test_SNR_significance_output_may_2018.mat
        roughness_vals=[0 2];
        D3b=D3b(:,:,[1 3]);
    end
end
clear Fgood F95 F98 Fbad
[sigma_good_98, iqr_good_98, iqr_scaled_good_98]=deal(NaN(numel(tau_values), numel(BGR_values), numel(roughness_vals)));
for kR=1:numel(roughness_vals)
    for kTau=1:numel(tau_values)
        for kB=1:numel(BGR_values)
            D3b{kTau, kB, kR}.SNR_significance=...
                interpn(SNR_F_table.SNR, SNR_F_table.W_surface_window_initial, SNR_F_table.BGR, SNR_F_table.P_NoiseOnly, ...
                max(-5, min(5, D3b{kTau, kB, kR}.SNR)), ...
                max(3, min(max(SNR_F_table.W_surface_window_initial), Htot))*ones(size(D3b{kTau, kB, kR}.SNR)), ...
                max(min(SNR_F_table.BGR), min(1e7, D3b{kTau, kB, kR}.BGR)));
            good=abs(D3b{kTau, kB, kR}.h_mean)<1 & isfinite(D3b{kTau, kB, kR}.h_mean)& isfinite(D3b{kTau, kB, kR}.h_mean)~=0 & abs(D3b{kTau, kB, kR}.dh_fit_dx) < 2/10  ;
            
            accepted=D3b{kTau, kB, kR}.SNR_significance < 0.02 & isfinite(D3b{kTau, kB, kR}.h_mean)& isfinite(D3b{kTau, kB, kR}.h_mean)~=0 & D3b{kTau, kB, kR}.sigma_h_mean< 1;% & D3b{kTau, kB, kR}.h_robust_spread < 1 ;%& abs(D3b{kTau, kB, kR}.dh_fit_dx) < 1/10;
            if sum(accepted) > 16 
                sigma_good_98(kTau, kB, kR)=std(D3b{kTau, kB, kR}.h_mean(accepted));
                iqr_good_98(kTau, kB, kR)=iqr(D3b{kTau, kB, kR}.h_med(accepted))/2;
                
                sigma_signal=sqrt((D3b{kTau, kB, kR}.dh_fit_dx.^2)*7.5.^2 + (0.7e-9*1.5e8).^2);
                N_noise=57*D3b{kTau, kB, kR}.w_surface_window_final/1.5e8*BGR_values(kB);
                N_sig=max(0, D3b{kTau, kB, kR}.N_final-N_noise);
                sigma_photon_est=sqrt((N_noise.*(D3b{kTau, kB, kR}.w_surface_window_final*.287).^2+N_sig.*sigma_signal.^2)./D3b{kTau, kB, kR}.N_final);
                sigma_per_photon=max(max(sigma_photon_est, D3b{kTau, kB, kR}.h_robust_spread), D3b{kTau, kB, kR}.h_rms);
                sigma_est=sigma_per_photon.*D3b{kTau, kB, kR}.sigma_h_mean;   
                
                iqr_scaled_good_98(kTau, kB, kR)=iqr(D3b{kTau, kB, kR}.h_mean(accepted)./sigma_est(accepted))/2;
            end
            Fgood(kTau, kB, kR)=mean(good);
            accepted95=D3b{kTau, kB, kR}.SNR_significance < 0.05 & D3b{kTau, kB, kR}.sigma_h_mean<1;
            accepted98=D3b{kTau, kB, kR}.SNR_significance < 0.02 & D3b{kTau, kB, kR}.sigma_h_mean<1;
            F95(kTau, kB, kR)=mean(accepted95 & good);
            F98(kTau, kB, kR)=mean(accepted98 & good);
            FgoodRejected95(kTau, kB, kR)=mean((~accepted95) & good);
            FgoodRejected98(kTau, kB, kR)=mean((~accepted98) & good);
            bad=(abs(D3b{kTau, kB, kR}.h_mean)>1 ) & isfinite(D3b{kTau, kB, kR}.h_mean) & isfinite(D3b{kTau, kB, kR}.h_mean)~=0;
            Fbad(kTau, kB, kR)=mean(bad);
            FBadRejected95(kTau, kB, kR)=mean((~accepted95) & bad);
            FBadRejected98(kTau, kB, kR)=mean((~accepted98) & bad);
            FBadAccepted95(kTau, kB, kR)=mean(accepted95 & bad);
            FBadAccepted98(kTau, kB, kR)=mean(accepted98 & bad);
        end
    end
end




figure(4); clf; set(gcf,'units','inches','position', [2 2 8.5 6],'defaultaxesfontsize', 12,'papersize', [9,6.5])
hax=cheek_by_jowl(2, 3, [0.15 0.15 .7, 0.7]);
yticks=[0.01    0.1   1   10]';
yticks_s=[1 2 4 8]';
for k=1:2
    axes(hax(k,1)); imagesc(BGR_values/1e6, tau_values, log10(sigma_good_98(:,:,k)), 'alphadata', double(isfinite(sigma_good_98(:,:,k))));   
    caxis(log10(range(yticks))); % hb=colorbar('north'); set(hb,'ytick', log10(yticks),'yticklabel', num2str((yticks)));   
    axes(hax(k,2)); imagesc(BGR_values/1e6, tau_values, log10(iqr_good_98(:,:,k)), 'alphadata', double(isfinite(iqr_good_98(:,:,k))));  
    caxis(log10(range(yticks))); %hb=colorbar('north'); set(hb,'ytick', log10(yticks),'yticklabel', num2str((yticks)));
    axes(hax(k,3)); imagesc(BGR_values/1e6, tau_values, log10(iqr_scaled_good_98(:,:,k)), 'alphadata', double(isfinite(iqr_scaled_good_98(:,:,k)))); 
    caxis(log10(range(yticks_s))); %hb=colorbar('north'); set(hb,'ytick', log10(yticks_s),'yticklabel', num2str((yticks)));
end
pos=get(hax(1, 1),'position');
hb=axes('position',  [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace(log10(min(yticks)), log10(max(yticks)), 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', log10(yticks),'xticklabel', num2str(yticks),'xaxislocation','top','ytick', []);
xlabel(hb, '\sigma, m');

pos=get(hax(1, 2),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace(log10(min(yticks)), log10(max(yticks)), 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', log10(yticks),'xticklabel', num2str(yticks),'xaxislocation','top','ytick', []);
xlabel(hb, '\sigma_R, m');

pos=get(hax(1, 3),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace(log10(min(yticks_s)), log10(max(yticks_s)), 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', log10(yticks_s),'xticklabel', num2str(yticks_s),'xaxislocation','top','ytick', []);
xlabel(hb, '\sigma_R/\sigma_{est}');
set(hax(:,3),'yaxislocation','right');
xlabel(hax(2,2),'Background rate, MHz')
set(hax(1,:),'xticklabel','');
yticks=[3 1 0.3 0.1 ]';
set(hax,'ytick', -log(yticks/3),'yticklabel', num2str(yticks, '%3.1f'));
set(hax(:,2),'yticklabel','');
ylabel(hax(1,1),{'Roughness=0 m', 'Signal ph./pulse'  });
ylabel(hax(2,1),{'Roughness=2 m', 'Signal ph./pulse'  });
colormap(jet*.8+.2); 
set(hax,'color', [.7 .7 .7]);
labels={'a','b','c','d','e','f'};
count=0;
clear ht;
for kr=1:2
    for kc=1:3
        count=count+1;
        axes(hax(kr, kc)); ht(count)=text(min(get(gca,'xlim'))+0.05, min(get(gca,'ylim')), labels{count});
    end
end
set(ht, 'HorizontalAlignment','left','VerticalAlignment','top','color','k','backgroundcolor','w','fontsize', 14);

%print -f4 -djpeg99 figures/linear_fit_accuracy

% map of the fraction of good segments rejected at each confidence level
figure(1); clf; set(gcf,'units','inches','position', [2 2 8.5 6],'defaultaxesfontsize', 12,'papersize', [9,6.5])
yticks=[0 0.1 0.2 0.3]'*100;
hax=cheek_by_jowl(2,3, [0.15 0.08 0.7 0.7]); 
for k=1:2
    axes(hax(k,1)); imagesc(BGR_values/1e6, tau_values, Fgood(:,:,k)); caxis([0.0 1]); 
    axes(hax(k,2)); imagesc(BGR_values/1e6, tau_values, FgoodRejected95(:,:,k)); caxis([0 0.3]);
    axes(hax(k,3)); imagesc(BGR_values/1e6, tau_values, FgoodRejected98(:,:,k)); caxis([0 0.3]);
end; figure(gcf); 

pos=get(hax(1, 1),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace(0, 100, 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', 0:25:100,'xticklabel', num2str((0:25:100)'),'xaxislocation','top','ytick', []);
xlabel(hb, '% valid');

pos=get(hax(1, 2),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace((min(yticks)), (max(yticks)), 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', (yticks),'xticklabel', num2str(yticks),'xaxislocation','top','ytick', []);
xlabel(hb, {'% valid, rejected' '@ 5% conf.'});

pos=get(hax(1, 3),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace((min(yticks)), (max(yticks)), 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', (yticks),'xticklabel', num2str(yticks),'xaxislocation','top','ytick', []);
xlabel(hb, {'% valid, rejected' '@2% conf.'});

set(hax(:,3),'yaxislocation','right');
xlabel(hax(2,2),'Background rate, MHz')
set(hax(1,:),'xticklabel','');
yticks=[3 1 0.3 0.1 ]';
set(hax,'ytick', -log(yticks/3),'yticklabel', num2str(yticks, '%3.1f'));
set(hax(:,2),'yticklabel','');
ylabel(hax(1,1),{'Roughness=0 m', 'Signal ph./pulse'  });
ylabel(hax(2,1),{'Roughness=2 m', 'Signal ph./pulse'  });

colormap(this_cmap)
figure(gcf);
%print -f1 -dpdf figures/f_rejected



% map of the fraction of bogus segments accepted at each confidence level,
% plus the f_bogus_accepted @ 2%.
figure(2); clf; set(gcf,'units','inches','position', [2 2 8.5 6],'defaultaxesfontsize', 12,'papersize', [9,6.5])
yticks=[0 0.1 0.2 0.3]'*100;
hax=cheek_by_jowl(2,4, [0.15 0.08 0.7 0.7]); 
for k=1:2
    axes(hax(k,1)); imagesc(BGR_values/1e6, tau_values, 100*Fbad(:,:,k)); caxis([0.0 100]); 
    axes(hax(k,2)); imagesc(BGR_values/1e6, tau_values, 100*FBadAccepted95(:,:,k)); caxis([0 30]);
    axes(hax(k,3)); imagesc(BGR_values/1e6, tau_values, 100*FBadAccepted98(:,:,k)); caxis([0 30]);
    axes(hax(k,4)); imagesc(BGR_values/1e6, tau_values, 100*FgoodRejected98(:,:,k)); caxis([0 30]);
end
 
pos=get(hax(1, 1),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace(0, 100, 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', 0:25:100,'xticklabel', num2str((0:25:100)'),'xaxislocation','top','ytick', []);
xlabel(hb, '% blunders');

pos=get(hax(1, 2),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace((min(yticks)), (max(yticks)), 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', (yticks),'xticklabel', num2str(yticks),'xaxislocation','top','ytick', []);
xlabel(hb, {'% blunders' '@ 5% conf.'});

pos=get(hax(1, 3),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace((min(yticks)), (max(yticks)), 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', (yticks),'xticklabel', num2str(yticks),'xaxislocation','top','ytick', []);
xlabel(hb, {'% blunders' '@2% conf.'});

pos=get(hax(1, 4),'position');
hb=axes('position',   [pos(1)+0.15*pos(3), pos(2)+1.05*pos(4), pos(3)*0.7 pos(4)*0.1]);
temp=linspace((min(yticks)), (max(yticks)), 250);
imagesc(temp, [0 1], temp); 
set(hb,'xtick', (yticks),'xticklabel', num2str(yticks),'xaxislocation','top','ytick', []);
xlabel(hb, {'% valid, rejected' '@2% conf.'});


set(hax(:,4),'yaxislocation','right');
xl=xlabel(hax(2,2),'Background rate, MHz');
set(xl,'position', get(xl,'position').*[2 1 1]);
set(hax(1,:),'xticklabel','');
set(hax(2,:),'xtick', 2:2:8);
yticks=[3 1 0.3 0.1 ]';
set(hax,'ytick', -log(yticks/3),'yticklabel', num2str(yticks, '%3.1f'));
set(hax(:,2:3),'yticklabel','');
ylabel(hax(1,1),{'Roughness=0 m', 'Signal ph./pulse'  });
ylabel(hax(2,1),{'Roughness=2 m', 'Signal ph./pulse'  });

colormap(my_rgb_cpt(128).*repmat(linspace(0.5, 1, 128)', [1, 3]))

labels={' a',' b',' c',' d',' e',' f',' g',' h'};
k=0;
for kr=1:2
    for kc=1:4 
        k=k+1;
        axes(hax(kr, kc)); ht(k)=text(min(get(gca,'xlim')), min(get(gca,'ylim')), labels{k});
    end
end
set(ht, 'HorizontalAlignment','left','VerticalAlignment','top','color','w','fontsize', 14,'fontweight','bold');


figure(gcf);
%print -f2 -dpdf figures/f_blunder_plus


figure(3);  clf
hax=cheek_by_jowl(3,3, [0.15 0.15 0.7 0.7]); 
for k=1:3
    axes(hax(k,1)); imagesc(BGR_values/1e6, tau_values, Fbad(:,:,k)); caxis([0. 1]); 
    axes(hax(k,2)); imagesc(BGR_values/1e6, tau_values, FBadRejected95(:,:,k)); caxis([0 1]);
    axes(hax(k,3)); imagesc(BGR_values/1e6, tau_values, FBadRejected98(:,:,k)); caxis([0 1]);
end; figure(gcf);

axes(hax(1,1)); title('F bad'); hb=colorbar('north');
axes(hax(1,2)); title('F bad & rejected @95%'); hb=colorbar('north');
axes(hax(1,3)); title('F bad & rejected @98%'); hb=colorbar('north');
axes(hax(2,1)); ylabel('optical thickness'); 
axes(hax(3,2)); xlabel('BGR, MHz');
colormap(this_cmap)
 
% INSIGHT: need to propagate the rms height, not the robust spread.  
