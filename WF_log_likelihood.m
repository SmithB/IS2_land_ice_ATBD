function L=WF_log_likelihood(dx, z, Hwin, dXshot, BGR, WF,  z0, slope, sigma, SigRatio)

% L=WF_log_likelihood(dx, z, Hwin, dXshot, BGR, WF, sigma, z0, slope, Nsignal)
c=299792458.; 

SigRatio=max(0, min(1, SigRatio));

% figure out how much of the power of the return is lost due to H
% windowing
dt_window_bottom= ((z0-dXshot*slope)-Hwin(1))/(c/2);
dt_window_top=((z0-dXshot*slope)-Hwin(2))/(c/2);

C_sig_top=interpn(WF.t,  WF.sigma, WF.C, dt_window_top, sigma*ones(size(dt_window_top))); 
C_sig_top(dt_window_top < min(WF.t))=0;
C_sig_top(dt_window_top > max(WF.t))=1;

C_sig_bottom=interpn(WF.t,  WF.sigma, WF.C,  dt_window_bottom, sigma*ones(size(dt_window_top))); 
C_sig_bottom(dt_window_bottom < min(WF.t))=0;
C_sig_bottom(dt_window_bottom > max(WF.t))=1;

F_sig=mean(C_sig_bottom-C_sig_top);

Ntot=length(dx);


% !!!!! REDO with SigRatio instead of SNR !!!
% assume: SNR=PPP Np / BHNp
% F = fraction of signal in window (=F_sig);
% Ns = BHNp F SNR
%    = Nn F SNR
% Ntot= Nn+ Nn F SNR
%     = Nn (1 + F SNR)
% Nn = Ntot/(1 + F SNR)
% Ns = Nn F SNR  = Nt F SNR / (1+ F SNR) = Nt / (1/(F SNR) + 1)
% SNR=large, F=1 -> Ns ~ Nt


Nnoise=Ntot./(1+F_sig.*SNR);
Nshots=length(dXshot);
Nsig=Ntot.*F_sig.*SNR./(1+F_sig.*SNR);


% Psig is the probability per unit time of seeing a photon within the
% window.  Integrated over the window, it should be equal to SNR
% Normalize the probabilities by 1/the window length, so that they're not 
% something ridiculous
Psig=Nsig/Ntot*interpn(WF.t, WF.sigma, WF.p, (z-z0-dx*slope)/(-c/2), max(WF.sigma(1),sigma))*diff(Hwin/(c/2));
Psig(~isfinite(Psig))=0;

% Pnoise is the probability of seeing a photon per unit time per shot.
% integrated over the window it should equal 1-SNR
% Pnoise=Nnoise/Ntot/diff(Hwin/(c/2)); % --- but we normalize by 1/ the window
% length
Pnoise=Nnoise/Ntot;


L=sum(log(Psig+Pnoise))+log(poisson_pdf(Nnoise, BGR*diff(Hwin)/(c/2)*Nshots));

% in the noise-only case---
% Pnoise=1/diff(Hwin/(c/2)) = 1/Dt
% L = Nnoise log(1/Dt) + log(poisson_pdf(Nnoise, Nnoise1))
%  This has a maximum at Nnoise1=Nnoise
%
%  log(poisson_pdf(k+delta,k))-log(poisson_pdf(k,k) ~ (log(k+delta)-log(k))*sqrt(delta)  
% so 
%  L1 (Nnoise+delta) = 
%  L+delta log(1/Dt) + sqrt(delta) (log((k+delta)/k)))


if false;
    dt=2.5e-9;
    histogram((z-z0-dx*slope)/(-c/2), -2e-7:dt:2e-7);
    hold on;
    plot(WF.t, Nsignal*N_sig_norm*dt*interpn(WF.t, WF.sigma, WF.p, WF.t, sigma)*dt + Pnoise*Nshots*dt,'r','linewidth', 3);
end
    

