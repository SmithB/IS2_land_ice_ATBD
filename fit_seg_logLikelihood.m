function [Zctr, mctr, sigctr, fctr, Lbest]=fit_seg_logLikelihood(D2, params, WF_data, BGR) 
 
x0_vals=(1:57)*0.7; x0_vals=x0_vals-20;
HW=params.H_window;
Zctr=0; mctr=0;
sigctr=mean(range(WF_data.sigma)); Fsig=0.5;
Nz=ceil(diff(range(D2.h))/(sigctr*1.5e8)*4);
ZR=linspace(min(D2.h), max(D2.h), Nz);
sigR=linspace(min(WF_data.sigma), max(WF_data.sigma), 20);
slopeR=linspace(-0.3, 0.3, 20);
fR=0:0.1:1;

% get the initial slope and z values with a grid search?
  
for k=1:5
    [Zctr, search_hist]=linesearch_q5(@(z0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, z0, mctr,sigctr,  Fsig), ZR,struct('xtol', 10.^-k));
    ZR=Zctr+[-sigctr sigctr]*1.5e8;
    [mctr, search_hist]=linesearch_q5(@(m0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, m0,sigctr,  Fsig), slopeR, struct('xtol', 0.1*10^-k));
    slopeR=mctr+[-1 1]*max([0.01, abs(mctr)/4]);
    [sigctr, search_hist]=linesearch_q5(@(sig0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr, sig0, Fsig), sigR, struct('xtol', 0.1*1e-9*10^-k));
    sigR=sigctr+[-1 1]*max(0.68e-9/2, sigctr*0.25);
    [fctr, search_hist]=linesearch_q5(@(Fsig0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, Zctr, mctr, sigctr,  Fsig0), fR, struct('xtol', max(0.01, 10^-k))) ;
    fR=min(1, [0.5 1.5]*fctr);
end

% do a last search for Z0
[Zctr, search_hist]=linesearch_q5(@(z0)-WF_log_likelihood(D2.x0-20, D2.h, [-1 1]*HW/2, x0_vals, BGR, WF_data, z0, 0,sigctr,  Fsig), ZR,struct('xtol', 0.001));
Lbest=max(-search_hist.z);


