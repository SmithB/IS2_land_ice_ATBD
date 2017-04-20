%function [m, sigma_m]=fit_seg_logLikelihood(x, z, x0, BGR, WF)

if ~exist('WF_data','var');
    load WF_est.mat
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
    
    for k=1:size(WF_data.p,2)
        WF_data.C(:,k)=cumtrapz(WF_data.t, WF_data.p(:,k));
    end
end

% make some test data:
mtrue=0.0125;
Ntrue=3*exp(0);
sigtrue=sqrt((mtrue*7.5/(1.5e8))^2+(0.68e-9)^2);
Ztrue=mtrue*20;
BGR=10e6;
HW=120;
ftrue=Ntrue*57/(Ntrue*57 + BGR*HW*57/(1.5e8));

%[D2, params]=make_ATL03_data(N_pulses, N_chan, roughness, N_per_pulse, WF, BGR, H_window, AT_slope)
[D2 , params]=make_ATL03_data(57, 4, 0, Ntrue, WF, BGR, HW+20, mtrue);
D2=index_struct(D2, abs(D2.h)<HW/2);
params.H_window=HW;

sigma_wide=0.25*params.sigma_x/1.5e8;
x0_vals=(1:57)*0.7; x0_vals=x0_vals-20;

%WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Ztrue, mtrue, sigtrue,  ftrue)

Zctr=0; mctr=0; sigctr=mean(range(WF_data.sigma)); Fsig=0.5;
figure(1); clf
Nnoise=BGR*40*57/(1.5e8);
Nz=ceil(diff(range(D2.h))/(sigctr*1.5e8)*4);
ZR=linspace(min(D2.h), max(D2.h), Nz);
sigR=linspace(min(WF_data.sigma), max(WF_data.sigma), 20);
slopeR=linspace(-0.3, 0.3, 20);
fR=0:0.1:1;
for k=1:20
    subplot(4,1,1);hold on;
    %fplot(@(z0)WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, z0, mctr, sigctr, Fsig), ZR([1 end])); xlabel('z0');
    [Zctr, search_hist]=linesearch_q5(@(z0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, z0, mctr, sigctr,  Fsig), ZR,struct('xtol', 0.001));
    plot(search_hist.x, -search_hist.z,'.');
    ZR=Zctr+[-sigctr sigctr]*1.5e8;
    plot(Zctr, WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr,sigctr,  Fsig),'*');

    subplot(4,1,2);hold on;
    %fplot(@(m0)WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, m0, sigctr,  Fsig), [-0.3 0.3]); xlabel('m');
    [mctr, search_hist]=linesearch_q5(@(m0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, m0,sigctr,  Fsig), slopeR, struct('xtol', 0.0001));
    plot(search_hist.x, -search_hist.z,'.');
    slopeR=mctr+[-1 1]*max([0.01, abs(mctr)/4]);
    plot(mctr, WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr,sigctr,  Fsig),'*');
      
    subplot(4,1,3);hold on;
    %fplot(@(sig0)WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr, sig0,  Fsig), [0.5 2]*sigctr); xlabel('sigma');
    [sigctr, search_hist]=linesearch_q5(@(sig0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr, sig0, Fsig), sigR, struct('xtol', 0.01*1e-9));
    plot(search_hist.x, -search_hist.z,'.');
    sigR=sigctr+[-1 1]*max(0.68e-9/2, sigctr*0.25);
    plot(sigctr, WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr,sigctr,  Fsig),'*');

    subplot(4,1,4);hold on;
    %fplot(@(Fsig0)WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr,sigctr,  Fsig0), [0 1]); xlabel('Fsig');
    [Fsig, search_hist]=linesearch_q5(@(Fsig0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr, sigctr,  Fsig0), fR, struct('xtol', 0.01)) ;
    plot(search_hist.x, -search_hist.z,'.');
    fR=min(1, [0.5 1.5]*Fsig);
    plot(Fsig, WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr, sigctr,  Fsig),'*');

end
subplot(4,1,1); axis tight; plot(Ztrue*[1 1], get(gca,'ylim'), 'k--');
subplot(4,1,2); axis tight; plot(mtrue*[1 1], get(gca,'ylim'), 'k--');
subplot(4,1,3); axis tight; plot(sigtrue*[1 1], get(gca,'ylim'), 'k--');
subplot(4,1,4); axis tight; plot(ftrue*[1 1], get(gca,'ylim'), 'k--');

figure;WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr, sigctr,  Fsig, true)




if false
    
    mtrue=0.0125;
    Ntrue=3*exp(-2);
    sigtrue=sqrt((mtrue*7.5/(1.5e8))^2+(0.68e-9)^2);
    Ztrue=mtrue*20;
    BGR=1e7;
    HW=120;
    ftrue=Ntrue*57/(Ntrue*57 + BGR*HW*57/(1.5e8));
    figure;
    clear Zctr1 mctr1 sigctr1 Nctr1 Lbest
    
    for k=1:100
        [D2 , params]=make_ATL03_data(57, 4, 0, Ntrue, WF, BGR, HW+20, mtrue);
        D2=index_struct(D2, abs(D2.h)<HW/2);
        params.H_window=HW;
        [Zctr1(k), mctr1(k), sigctr1(k), Fctr1(k), Lbest(k)]=fit_seg_logLikelihood(D2, params, WF_data, BGR);
        if abs(Zctr1(k)-Ztrue)> 1; break; end
    end
    subplot(2,2,1); histogram(Zctr1, 20); xlabel('z0');
    subplot(2,2,2); histogram(sigctr1); xlabel('sigma');
    subplot(2,2,3); histogram(mctr1); xlabel('m');
    subplot(2,2,4); histogram(Fctr1); xlabel('N');
end
 