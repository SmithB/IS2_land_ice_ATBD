BGR=1e7;
sigma_rx=0.68e-9:0.1e-9:(1.5/1.5e8);
%Nsigma=2:.1:10;
Nsigma=3;
clear sigma_med sigma_bar
%sigma_rx=0.5/1.5e8
kN=1;
%for kN=1:length(Nsigma)

for k=1:length(sigma_rx)
    HW=max(3,2*Nsigma(kN)*sigma_rx(k)*1.5e8);
    Nnoise=ceil(HW*BGR/1.5e8*57);
    h=[(rand(Nnoise*2, 10240)-0.5)*HW*2; sigma_rx(k)*1.5e8*randn(57*3, 10240)];
    Hctr=h-mean(h);
    h(Hctr >  HW/2)=NaN;
    h(Hctr < -HW/2)=NaN;
    
    
    sigma_med(k, kN )=std(nanmedian(h));
    sigma_bar(k, kN)=std(nanmean(h));
end

    
clf; plot(sigma_med); hold on; plot(sigma_bar); legend('median','mean')