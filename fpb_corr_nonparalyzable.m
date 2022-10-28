function [med, centroid, count, N_fpb_corr, sigma_med, sigma_centroid, minGain, gain_est]=...
    fpb_corr_nonparalyzable(t_WF_full, counts,  N_chan, N_pulses, t_dead, dt)
c=3e8;
t_WF_full=t_WF_full(:);

if ~exist('signal_threshold_for_gain_corr','var')
    signal_threshold_for_gain_corr=0.02;  % PE/dead time/channel/pulse.  At 10 MHz noise, get .00267, so threshold of .01 is between 2 and 3 sigma
end

N0_full=counts/N_chan/N_pulses;

% make sure nothing funny happens at the start and end of the WF
N0_full(1)=0; 
N0_full(end)=0;
N_dt_bins=floor(t_dead/dt);  


% for each detected photon, the pixel is inactive for the next N_dt
% bins. Since there are N_pulses*N_chan detector channels in the array, the
% resulting gain is:
K=[zeros(N_dt_bins,1); ones(N_dt_bins,1)];
gain_est=(1-conv(counts(:), K(:), 'same')/(N_pulses*N_chan));

% disp('WARNING: running ASAS check!!!');
% % test vs. ASAS code:
% N_dt_ASAS=ceil(t_dead/dt);
% p_dead=zeros(size(counts));
% N_PC=N_pulses*N_chan;
% for k=2:length(p_dead)
%     ii=max(1, k-N_dt_ASAS):k-1;
%     p_dead(k)=sum(counts(ii))/N_PC;
% end
% gain_ASAS=1-p_dead;


minGain=min(gain_est);
 
% renormalize to histogram counts
N_fpb_corr=counts(:)./gain_est(:);

if minGain < 0.1
     [med, centroid, count, sigma_med, sigma_centroid]=deal(NaN);
     return
end

% calculate the stats AFN height
[med, centroid, count, sigma_med, sigma_centroid]=calc_stats(counts(end:-1:1), gain_est(end:-1:1), t_WF_full(end:-1:1)*(-1.5e8) );

%----------------------------------------------------
function [med, centroid, N, sigma_med, sigma_centroid]=calc_stats(WF, gain, t_WF)

if all(WF==0) 
    [med, centroid, N, sigma_med, sigma_centroid]=deal(NaN);
    return
end

WFc=WF(:)./gain(:);
centroid=sum(t_WF(:).*WFc(:))./sum(WFc(:));
N=sum(WF./gain);

t_40_50_60=percentile_of_histogram([0.4 0.5 0.6], t_WF, WF./gain);
med=t_40_50_60(2);

sigma_WF=sqrt(WF)./gain;
sigma_CDF_med=sqrt([0.5*sum(sigma_WF.^2)./N^2]);
sigma_med=abs(diff(t_40_50_60([1 3])))/0.2*sigma_CDF_med;

sigma_centroid=sqrt(sum((WF./gain.*(t_WF(:)-centroid)/N).^2));


%----------------------------------------------------
function X=percentile_of_histogram(P, bins, counts)
bin_width=[diff(bins(:)); bins(end)-bins(end-1)];
edges=[bins(1)-bin_width(1)/2; bins(:)+bin_width/2];
C=[0; cumsum(counts(:))]; C=C/C(end);

for k=1:length(P)
    % check if any C are equal to P(k)
    eq_els=abs(C-P(k)) < 100*eps;
    if any(eq_els)
        X(k)=mean(edges(eq_els));
        continue
    end
    i_minus=find(C<P(k)-100*eps, 1, 'last');
    i_plus=find(C>=P(k)+100*eps, 1, 'first');   
    if C(i_plus)==C(i_minus)
        X(k)=0.5*(edges(i_plus)+edges(i_minus));
    else
        X(k)=((P(k)-C(i_minus))*edges(i_plus)+(C(i_plus)-P(k))*edges(i_minus))/(C(i_plus)-C(i_minus));
    end
end

