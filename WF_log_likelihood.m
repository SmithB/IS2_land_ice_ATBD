function [L, eL, sigmaL]=WF_log_likelihood(dx, z, Hwin, dXshot, BGR, WF,  z0, slope, sigma, SigRatio, DoPlot)

% vectorize.  (not really- just build a loop)
[z01, slope1, sigma1, SigRatio1]=ndgrid(z0(:), slope(:), sigma(:), SigRatio(:));
z01=z01(:); slope1=slope1(:); sigma=sigma1(:); SigRatio=SigRatio1(:);

% L=WF_log_likelihood(dx, z, Hwin, dXshot, BGR, WF, sigma, z0, slope, Nsignal)
c=299792458.;
Twin=diff(Hwin)/(c/2);
 
% figure out how much of the power of the return is lost due to H
% windowing
for k=1:length(slope1(:))
    
    [z0, slope, sigma, SigRatio]=deal(z01(k), slope1(k), sigma1(k), SigRatio1(k));
    if false
        dt_window_bottom= ((z0-dXshot*slope)-Hwin(1))/(c/2);
        dt_window_top=((z0-dXshot*slope)-Hwin(2))/(c/2);
        
        C_sig_top=interpn(WF.t,  WF.sigma, WF.C, dt_window_top, sigma*ones(size(dt_window_top)));
        C_sig_top(dt_window_top < min(WF.t))=0;
        C_sig_top(dt_window_top > max(WF.t))=1;
        
        C_sig_bottom=interpn(WF.t,  WF.sigma, WF.C,  dt_window_bottom, sigma*ones(size(dt_window_top)));
        C_sig_bottom(dt_window_bottom < min(WF.t))=0;
        C_sig_bottom(dt_window_bottom > max(WF.t))=1;
        
        Fsig=mean(C_sig_bottom-C_sig_top);
    end
    Fsig=1;
    Ntot=length(dx);
    
    Nsig=Fsig.*SigRatio.*Ntot;
    Nnoise=Ntot-Nsig;
    Nshots=length(dXshot);
    
    % Psig is the probability per unit time of seeing a photon within the
    % window.  Integrated over the window, it should be equal Fsig
    % Normalize the probabilities by 1/the window length, so that they're not
    % something ridiculous
    Psig=Nsig.*interpn(WF.t, WF.sigma, WF.p, (z-z0-dx*slope)/(-c/2), max(WF.sigma(1),sigma))*Twin;
    Psig(~isfinite(Psig))=0;
    
    % Pnoise is the probability of seeing a photon per unit time per shot.
    % integrated over the window it should equal 1-SNR
    % Pnoise=Nnoise*diff(Hwin/(c/2)); % --- but we normalize by 1/ the window
    % length
    Pnoise=Nnoise; 
    L(k)=sum(log(Psig+Pnoise))+log(poisson_pdf(Nnoise, BGR*diff(Hwin)/(c/2)*Nshots));
    
    % in the noise-only case---
    % Pnoise=1/diff(Hwin/(c/2)) = 1/Dt
    % L = Nnoise log(1/Dt) + log(poisson_pdf(Nnoise, Nnoise1))
    %  This has a maximum at Nnoise1=Nnoise
    %
    %  log(poisson_pdf(k+delta,k))-log(poisson_pdf(k,k) ~ (log(k+delta)-log(k))*sqrt(delta)
    % so
    %  L1 (Nnoise+delta) =
    %  L+delta log(1/Dt) + sqrt(delta) (log((k+delta)/k)))
    
    
    if exist('DoPlot','var') && DoPlot
        clf;
        dt=3e-9;
        Twin=diff(Hwin)/(c/2);
        TT=-Twin/2:dt:Twin/2;

        histogram((z-z0-dx*slope)/(-c/2), TT);
        hold on;
        P1=interpn(WF.t, WF.sigma, WF.p, TT, sigma); P1(~isfinite(P1))=0;
        plot(TT, Nsig*P1*dt + Nnoise*dt/Twin,'r','linewidth', 3);
        figure(gcf);
    end
end


% Note that the expected value for LL (i.e. the best-fitting model)
% for a random distribution of photons
% over the window is N*mean(LL)+poisson_PDF(N, N);
% the standard deviation is std(LL) (right?)
% can check whether L_recovered > E(L)+M*sigma



% What is the sampling standard deviation for LL for a true model? (LISN)
 