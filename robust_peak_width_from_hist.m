function [sigma_hat, med]=robust_peak_width_from_hist(x, z, N_BG, XR)

% C_tot=C_sig+C_noise
%      =C_sig(t)+(t-t0)BGR
% want t1 such that C_sig(t1)=f
% Ctot(t1)=f+(t1-t0)BGR
% look for first point where Ctot > 0.16 Nsig + (t1-t0)*BGR
% and last point where Ctot < 0.84 Nsig + (t1-t0)*BGR
% the difference gives the peak width

if ~exist('TR','var')
    XR=[min(x) max(x)];
end
BGR=N_BG/diff(XR);

% what if this doesn't work? -- if N_sig < sqrt(N_BG) this will fail about
% half the time.  
N_sig=sum(z)-N_BG;

C=cumsum(z);

i0=find(C<0.16*N_sig + (x-XR(1))*BGR, 1, 'last');
if isempty(i0); i0=1;end
i1=find(C>0.84*N_sig + (x-XR(1))*BGR, 1, 'first');
if isempty(i1); i1=length(x);end
sigma_hat=abs(diff(x([i0 i1]))/2);

if nargout==2;
    i0=find(C<0.5*N_sig + (x-XR(1))*BGR, 1, 'last');
    if isempty(i0); i0=1;end
    i1=find(C>0.5*N_sig + (x-XR(1))*BGR, 1, 'first');
    if isempty(i1); i1=length(x);end
    if i0~=i1
        med=interp1(C([i0 i1])-(x([i0 i1])-XR(1))*BGR, x([i0 i1]), 0.5*N_sig);
    else
        med=x(i0);
    end
end

% test:
if false
    clear
    load ATLxx_example/WF_est;
    sigma_0=0.5*diff(wf_percentile(WF.t, WF.p, [0.16 0.84]));
    % the minimum SNR is 10 PE /( 6 MHz * 20 m /1.5e8 m/s * 57 pulses)=.13
    
    % need a quick test.  Do this for
    %-- signal between 10 and 684 PE/shot
    %-- BG between 0.1 and 10 MHz;
    %-- sigma between sigma_0 and 10*sigma_0;
    
    
    % set Nrate to 1e7 * 57
    N_pulses=57;
    Nrate=N_pulses*8e6;
    N_MC=1024;
    H_window=20;
    
    
    dz=2.5e-3;
    t_hist=(-H_window:dz:H_window/2)' /1.5e8;
    sigma_vals=sigma_0*(1+logspace(-2, 1, 10));
    Nsignal_vals=logspace(log10(10), log10(700), 21);
    N_noise=H_window/1.5e8*Nrate;
    [h.sigma.bar, h.sigma.sigma_hat, h.sigma.E95, h.med.bar, h.med.sigma_hat, h.med.E95]=deal(NaN(length(Nsignal_vals), length(sigma_vals)));
    for kNS=1:length(Nsignal_vals)   % LOOP OVER SIGNAL COUNT
        for kSig=1:length(sigma_vals)  % LOOP OVER PEAK WIDTH
            sigma_extra=sqrt(sigma_vals(kSig).^2-sigma_0.^2);
            [sigma_hat, med]=deal(NaN(N_MC,1));
            for kMC=1:N_MC  % MONTE CARLO LOOP
                % set the signal strength to Nrate*SNR*(HW/1.5e8)
                this_N_sig=Nsignal_vals(kNS);
                t=random_wf_sim(WF.t, WF.p, poisson_rv(this_N_sig, 1));
                t=t+sigma_extra*randn(size(t));
                t=sort([t(:); H_window/1.5e8*2*(rand(poisson_rv(N_noise,1),1)-0.5)]);
                
                
                P_hist=my_histc(t, t_hist);
                this_sigma_hat=H_window/1.5e8/6;
                  %for kk=1:5
                
                signal_PE_count_est=length(t)-N_noise;
                if signal_PE_count_est < sqrt(length(t));  % trouble if the signal level is poor
                    bins_2p5=2.5/dz;
                    
                    P1=conv(P_hist(:), ones(bins_2p5,1),'same');
                    N_noise_expected=bins_2p5*dz*Nrate;
                    els=P1>2*N_noise_expected;
                    els(min(find(els)):max(find(els)))=true;
                else
                    els=true(size(t_hist));
                end
                
                [this_sigma_hat, this_med]=robust_peak_width_from_hist(t_hist(els), P_hist(els), N_noise, 6*this_sigma_hat);
                els=abs(t_hist-this_med) < 3*this_sigma_hat;
                %    if kk==1; sigma_hat_initial(kMC)=this_sigma_hat; end
                %end
                sigma_hat(kMC)=this_sigma_hat;
                med(kMC)=this_med;
            end
            h.sigma.bar(kNS, kSig)=mean(sigma_hat);
            h.sigma.sigma_hat(kNS, kSig)=iqr(sigma_hat)/2;
            h.sigma.E95(kNS, kSig)=diff(percentile(sigma_hat, [0.025 0.975]))/2;
            h.med.bar(kNS, kSig)=mean(med);
            h.med.sigma_hat(kNS, kSig)=iqr(med)/2;
            h.med.E95(kNS, kSig)=diff(percentile(med, [0.025 0.975]))/2;
        end
    end
    
    % display:  Errors in sigma, ranges in sigma
    figure;
    clf;
    subplot(1,3,1);
    surf(sigma_vals*1.5e8,Nsignal_vals,  1.5e8*(h.sigma.bar-repmat(sigma_vals(:)',[length(Nsignal_vals), 1]))); colorbar; xlabel('sigma, m'); ylabel('SNR');
    set(gca,'yscale','log'); xlabel('pulse width, m'); ylabel('signal PE'); zlabel('bias in return width');
    subplot(1,3,2);
    surf(sigma_vals*1.5e8, Nsignal_vals, h.sigma.sigma_hat*1.5e8); colorbar;
    set(gca,'yscale','log'); xlabel('pulse width, m'); ylabel('signal PE'); zlabel('68% range in return width');
    
    subplot(1,3,3);
    surf( sigma_vals*1.5e8, Nsignal_vals, h.sigma.E95*1.5e8); colorbar;set(gca,'yscale','log');
    set(gca,'yscale','log'); xlabel('pulse width, m'); ylabel('signal PE'); zlabel('95% range in return width');
    
    h.sigma_vals=sigma_vals; h.Nsignal_vals=Nsignal_vals;
    
    save('rpw_test_8Mhz', 'h');
end
