function [med, centroid, count, N_fpb_corr, sigma_med, sigma_centroid, minGain, gain_full]=fpb_corr(t_WF_full, counts,  N_chan, N_pulses, t_dead, signal_threshold_for_gain_corr, dt)
c=3e8;

if ~exist('signal_threshold_for_gain_corr','var')
    signal_threshold_for_gain_corr=0.02;  % PE/dead time/channel/pulse.  At 10 MHz noise, get .00267, so threshold of .01 is between 2 and 3 sigma
end

%dt=diff(t_WF_full(1:2));

N0_full=counts/N_chan/N_pulses;

% make sure nothing funny happens at the start and end of the WF
N0_full(1)=0; 
N0_full(end)=0;
N_dt_bins=floor(t_dead/dt);  % should be -1.  To match ASAS, I need +3.  The 'floor' is a problem here- sometimes rounds down by a full bin.  

% calculate the number of Ph per dead time
N_per_dt=conv(N0_full,  [zeros(N_dt_bins,1); ones(N_dt_bins,1)],'same');

% only calculate the gain for bins for which the PH/deadtime rate is greater than 0.1 
if ~any(N_per_dt>signal_threshold_for_gain_corr)
    % generate a default gain estimate for the full waveform
    % If we haven't calculated a gain value, assume it's equal to 1t
    gain_full=ones(size(N0_full));
     
    N_fpb_corr=N0_full./gain_full*N_pulses*N_chan;
    [med, centroid, count, sigma_med, sigma_centroid]=calc_stats(N0_full(:)*N_pulses*N_chan, gain_full(:), t_WF_full(:)*(-1.5e8) );
    minGain=1;
    return
end
TR=range(t_WF_full(N_per_dt>signal_threshold_for_gain_corr))+[-1 1]*t_dead;
gain_calc_bins=t_WF_full >=TR(1) & t_WF_full <= TR(2);
N0=N0_full(gain_calc_bins);
 
%N_dt_bins=floor(t_dead/dt)-1;

% initialize the correction
N=N0;
N_corr=N0;
iteration=0;
last_gain=ones(size(N_corr));
gain=zeros(size(N_corr));

% iterate the solution until the gain converges
while max(abs(last_gain-gain)) > 0.001 && iteration < 20
    iteration=iteration+1;
    last_gain=gain;
    %for iteration=1:10;
    gain=exp(conv(log(max(10000*eps, 1-N_corr)), [zeros(N_dt_bins,1); ones(N_dt_bins, 1)],'same'));
    N_corr=N./gain;
end
     
% generate a default gain estimate for the full waveform
% If we haven't calculated a gain value, assume it's equal to 1
gain_full=ones(size(N0_full));
gain_full(gain_calc_bins)=gain;
minGain=min(gain_full);

% renormalize to histogram counts
N_fpb_corr=N0_full./gain_full*N_pulses*N_chan;
N0_full=N0_full*N_pulses*N_chan;
if minGain < 0.1
     [med, centroid, count, sigma_med, sigma_centroid]=deal(NaN);
     return
end

[med, centroid, count, sigma_med, sigma_centroid]=calc_stats(N0_full(:), gain_full(:), t_WF_full(:)*(-1.5e8) );

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

sigma_centroid=sqrt(sum((WF./gain.*(t_WF-centroid)/N).^2));


%----------------------------------------------------
function X=percentile_of_histogram(P, bins, counts)
bin_width=[diff(bins(:)); bins(end)-bins(end-1)];
edges=[bins(1)-bin_width(1)/2; bins(:)+bin_width/2];
C=[0; cumsum(counts(:))]; C=C/C(end);
for k=1:length(P)
    % check if any C are equal to P(k)
    eq_els=abs(C-P(k)) < 100*eps;
    if any(eq_els);
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

