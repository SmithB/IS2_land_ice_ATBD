function varargout=test_fpb_corr(varargin)

if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end


%------------------------------------------------------
function [D3u, R, HW]=make_test_data 
% run a tx-pulse-shape correction simulation for a range of surface roughnesses (R).
% output values give the median and mean height offset.
%
% inputs: 
%    WF:   The transmit pulse.  struct: t: time, p: power.  t should be offset so that the waveform centroid is zero
%
% outputs:
%   R: the roughness values for which the simulation was calculated, m
%   N_per_pulse: the signal strengths (photons per pulse) for the simulation
%   bar: means of output values, struct with fields:
%      med: corrected median photon height
%      centroid: corrected centroid heights
%      med_uncorr: uncorrected median photon heights
%      centroid_uncorr: uncorrected centroid photon heights
%   med: medians of output values, struct with same fields as bar

load WF_est;
% use a moderately-unsaturated strong WF.     
N_chan=12;
N_per_pulse=N_chan*1.25;
R=0:.125:1;
HW=[3:.5:8];
[R, HW]=meshgrid(R, HW);
% aligh the transmit waveform with time zero at its centroid.
dt_WF=(WF.t(2)-WF.t(1));
WF.t=WF.t-sum(WF.t.*WF.p)./sum(WF.p);
BGR=5e6; 

% to get an SNR of 100

load TX_shape_corr_table

% N.B.  This can be run as a parfor if you have the parallel tool box
parfor kk=1:numel(R)     
        fprintf(1, 'R=%d, HW=%3.1d\n', R(kk), HW(kk));
        
        % run this for enough iterations to beat down the noise in the means and medians caused by the spread of the RX pulse
        N_pulses=floor(5e3*(R(kk).^2+.24^2)/(.24^2));
        % Generate the fake data
        [D2, params]=test_tx_shape_corr('make_data',N_pulses, N_chan, R(kk), N_per_pulse, WF, BGR);
        
        % generate the uncorrected data;
        [D3u(kk)]=calc_h(D2, HW(kk), params);
end
D3u=reshape(D3u, size(R));

%====================
function [D3c, med_u, bar_u, med_c, bar_c]=correct_data(D3u, LUT)

[med_u, bar_u, med_c, bar_c]=deal(NaN(size(D3u)));
D3c=D3u;
for k=1:numel(D3u)
    dMed=interpn(LUT.W_pulse, LUT.HW, LUT.SNR, LUT.delta_med, ...
        min(max(D3u(k).sigma_rx*1.5e8, min(LUT.W_pulse)), max(LUT.W_pulse)), ...
        min(max(D3u(k).HW, min(LUT.HW)), max(LUT.HW)), ...
        min(max(D3u(k).SNR,  min(LUT.SNR)), max(LUT.SNR)));
    D3c(k).dMed=dMed;
    D3c(k).med=D3u(k).med+dMed;
    
    dMean=interpn(LUT.W_pulse, LUT.HW, LUT.SNR, LUT.delta_mean, ...
        min(max(D3u(k).sigma_rx*1.5e8, min(LUT.W_pulse)), max(LUT.W_pulse)), ...
        min(max(D3u(k).HW, min(LUT.HW)), max(LUT.HW)), ...
        min(max(D3u(k).SNR,  min(LUT.SNR)), max(LUT.SNR)));
    D3c(k).d_centroid=dMean;
    D3c(k).centroid=D3u(k).centroid+dMean;
    
    med_u(k)=mean(D3u(k).med);
    bar_u(k)=mean(D3u(k).centroid);
    med_c(k)=mean(D3c(k).med);
    bar_c(k)=mean(D3c(k).centroid);
    
end




%-------------------------------------------------------
function [D3u]=calc_h(D2, HW, params)

% calculate h from received-pulse data.  
% inputs:
% D2: structure array giving PE heights and pulse numbers
% HW: surface-window height, in meters, used to truncate the data
% params: a structured array that includes (at minimum) the transmit-waveform shape and the background rate
%
% outputs:
% D3: corrected heights,based on the truncated PE distribution,  including the TX-shape correction
% D3u: Uncorrected heights based on the truncated PE distribution not including the TX-shape correction
% D3a: heights based on all of the PE, with no correction and no truncation.

segs=1:58:max(D2.pulse);
[D3u.med, D3u.centroid, D3u.count, D3u.sigma_rx, D3u.SNR, D3u.HW]=deal(NaN(length(segs)-1,1));
 %D3u=D3;
CONSTANTS.dt_hist=25e-12;
if length(segs)==1
    segs(2)=segs(1)+58;
end
[t_bin, bin_count]=deal(cell(length(segs)-1,1));
for k=1:length(segs)-1
    in_bin=D2.pulse>=segs(k) & D2.pulse <=segs(k+1);
     
    % take an initial subset of pulses in the segment
    D2sub0=index_struct(D2, in_bin);
    
    % iterate 5x to find the truncated mean:
    Hmean=median(D2sub0.h);
    in_HW=abs(D2sub0.h-Hmean)<HW/2;
    for k_it=1:50
        Hmean=mean(D2sub0.h(in_HW));
        in_HW=abs(D2sub0.h-Hmean)<HW/2;
    end
    D2sub=index_struct(D2sub0, in_HW);
    
    % calc the uncorrected statistics;
    D3u.med(k)=median(D2sub.h);
    D3u.centroid(k)=mean(D2sub.h);
    
    D3u.sigma_rx(k)=robust_peak_width_CDF(D2sub.h/(-1.5e8), 58*HW/1.5e8*params.NoiseRate, HW*[-0.5 0.5]/(1.5e8));
    
    N_noise=58*HW/1.5e8*params.NoiseRate;
    N_signal_est=sum(in_HW)-N_noise;
    D3u.SNR(k)=N_signal_est/N_noise;
    D3u.HW(k)=HW;
    
end


