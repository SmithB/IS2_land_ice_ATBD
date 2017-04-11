%function [m, sigma_m]=fit_seg_logLikelihood(x, z, x0, BGR, WF)


load /Volumes/ice1/ben/sdt/WF_est.mat
t0=[-64:.1:64]'*1e-9;
WF.p=interp1(WF.t, WF.p, t0); 
WF.p(~isfinite(WF.p))=0;
WF.t=t0;
WF.t=WF.t-trapz(WF.t, WF.t.*WF.p)./trapz(WF.t, WF.p);


WF_data.sigma=(0.68:0.05:20)*1e-9; % ns
sigma_extra_vals=sqrt(WF_data.sigma.^2-(0.68e-9)^2); 
WF_data.p=NaN(length(WF.t), length(sigma_extra_vals));
WF_data.t=WF.t(:);

dt=diff(WF.t(1:2));
for k=1:length(sigma_extra_vals);
    Nsamps=sigma_extra_vals(k)/dt;
    if round(Nsamps) > 0;
        t_gaussian=[-ceil(Nsamps*3):ceil(Nsamps*3)]'*dt;
        WF_data.p(:, k)=conv(WF.p, dt*gaussian(t_gaussian, 0, sigma_extra_vals(k)),'same');
    else
        WF_data.p(:,k)=WF.p;
    end
    WF_data.p(:,k)=WF_data.p(:,k)./trapz(WF_data.t, WF_data.p(:,k));
end

for k=1:size(WF_data.p,2);
    WF_data.C(:,k)=cumtrapz(WF_data.t, WF_data.p(:,k));
end


% make some test data:
mtrue=0.25;
Ntrue=3;
sigtrue=sqrt((mtrue*7.5/(1.5e8))^2+(0.68e-9)^2);
Ztrue=mtrue*20;
%[D2, params]=make_ATL03_data(N_pulses, N_chan, roughness, N_per_pulse, WF, BGR, H_window, AT_slope)
[D2 , params]=make_ATL03_data(57, 4, 0, Ntrue, WF, BGR, 40, mtrue);

sigma_wide=0.25*params.sigma_x/1.5e8;
x0_vals=(1:57)*0.7; x0_vals=x0_vals-20;


Zctr=0; mctr=0; sigctr=mean(range(WF_data.sigma)); Nctr=1;
figure(1); clf
Nnoise=BGR*40*57/(1.5e8);
for k=1:3;
    subplot(4,1,1);hold on;
    fplot(@(z0)WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, z0, mctr, sigctr,  Nctr*57/Nnoise), [-10, 10]); xlabel('z0');
    Zctr=fminsearch(@(z0)-WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, z0, 0,sigma_wide,  Nctr*57/Nnoise), Zctr)
    plot(Zctr, WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, mctr,sigctr,  Nctr*57/Nnoise),'x');

    subplot(4,1,2);hold on;
    fplot(@(m0)WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, m0, sigctr,  Nctr*57/Nnoise), [-0.3 0.3]); xlabel('m');
    mctr=fminsearch(@(m0)-WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, m0,sigctr,  Nctr*57/Nnoise), mctr)
    plot(mctr, WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, mctr,sigctr,  Nctr*57/Nnoise),'x');
      
    subplot(4,1,3);hold on;
    fplot(@(sig0)WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, mctr, sig0,  Nctr*57/Nnoise), range(WF_data.sigma)); xlabel('sigma');
    sigctr=abs(fminsearch(@(sig0)-WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, mctr, sig0,  Nctr*57/Nnoise), sigctr))
    plot(sigctr, WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, mctr,sigctr,  Nctr*57/Nnoise),'x');

    subplot(4,1,4);hold on;
    fplot(@(N0)WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, mctr,sigctr,  N0*57/Nnoise), [0.5 4]); xlabel('N');
    Nctr=fminsearch(@(N0)-WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, mctr, sigctr,  N0*57/Nnoise), Nctr) 
    plot(Nctr, WF_log_likelihood(D2.x0-20, D2.h, [-20 20], x0_vals, BGR, WF_data, Zctr, mctr, sigctr,  Nctr*57/Nnoise),'x');

end
subplot(4,1,1); plot(Ztrue*[1 1], get(gca,'ylim'), 'k--');
subplot(4,1,2); plot(mtrue*[1 1], get(gca,'ylim'), 'k--');
subplot(4,1,3); plot(sigtrue*[1 1], get(gca,'ylim'), 'k--');
subplot(4,1,4); plot(Ntrue*[1 1], get(gca,'ylim'), 'k--');

figure; 
clear Zctr1 mctr1 sigctr1 Nctr1
for k=1:100; 
    [D2 , params ]=make_ATL03_data(57, 4, 0, 1, WF, 1e7, 40, .25); 
    [Zctr1(k), mctr1(k), sigctr1(k), Nctr1(k)]=fit_seg_logLikelihood(D2, params, WF_data, 1e7); 
end
subplot(2,2,1); histogram(Zctr1); xlabel('z0'); 
subplot(2,2,2); histogram(sigctr1); xlabel('sigma'); 
subplot(2,2,3); histogram(mctr1); xlabel('m'); 
subplot(2,2,4); histogram(Nctr1); xlabel('N');

 