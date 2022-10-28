function varargout=test_tx_shape_corr_bnds(varargin)

if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end

%------------------------------------------------------
function [R, HW, bar_out, med_out]=run_sim(WF)
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

% use a moderately-unsaturated strong WF.     
N_chan=12;
N_per_pulse=N_chan*1.25;
R=0:.125:1;
HW=[3:.5:8];
[R, HW]=meshgrid(R, HW);
% aligh the transmit waveform with time zero at its centroid.
dt_WF=(WF.t(2)-WF.t(1));
WF.t=WF.t-sum(WF.t.*WF.p)./sum(WF.p);
BGR=1e6; 
% N.B.  This can be run as a parfor if you have the parallel tool box
for kk=1:numel(R)     
        fprintf(1, 'R=%d, HW=%3.1d\n', R(kk), HW(kk));
        
        % run this for enough iterations to beat down the noise in the means and medians caused by the spread of the RX pulse
        N_pulses=floor(5e3*(R(kk).^2+.24^2)/(.24^2));
        % Generate the fake data
        [D2, params]=make_data(N_pulses, N_chan, R(kk), N_per_pulse, WF, BGR);
        
        % test the correction
        [D3, D3u, D3a]=calc_h(D2, HW(kk), params);
        % save output : calculate the median of each output parameter
        temp_med=median(D3.med);
        temp_centroid=median(D3.centroid);
        temp_med_uncorr=median(D3u.med);
        temp_centroid_uncorr=median(D3u.centroid);
        med(kk)=struct('med', temp_med,'centroid',temp_centroid,'centroid_uncorr', temp_centroid_uncorr, 'med_uncorr', temp_med_uncorr);
       
        % save output : calculate the mean of each output parameter
        temp_med=mean(D3.med);
        temp_centroid=mean(D3.centroid);
        temp_med_uncorr=mean(D3u.med);
        temp_centroid_uncorr=mean(D3u.centroid);
        bar(kk)=struct('med', temp_med,'centroid',temp_centroid,'centroid_uncorr', temp_centroid_uncorr, 'med_uncorr', temp_med_uncorr);
%    end
%end
end
% collect the output:

f=fieldnames(bar);
for kf=1:length(f); 
    med_out.(f{kf})=reshape([med.(f{kf})], size(R));
    bar_out.(f{kf})=reshape([bar.(f{kf})], size(R));
end

%-------------------------------------------------------
function out_dir=generate_test_case_data
load WF_est.mat
R=[0 .125 0.25 .5 .75 1 1.25 ];
HW=max(3, 2*R);
N_chan=12;
N_per_pulse=N_chan*1.25;
% aligh the transmit waveform with time zero at its centroid.
dt_WF=(WF.t(2)-WF.t(1));
WF.t=WF.t-sum(WF.t.*WF.p)./sum(WF.p);
out_dir=['/home/ben/Dropbox/projects/IS2_ATBD/tx_corr_test/data_' strrep(datestr(now),' ','-')];
if ~exist(out_dir,'dir'); mkdir(out_dir);end
ph_fields={'t_ph', 'h', 'pulse','channel'};
BGR=1e7;
N_pulses=57;
for count=1:1
    for kk=1:length(R);
        [D2, params]=make_data(N_pulses, N_chan, R(kk), N_per_pulse, WF, BGR);
        %D2=index_struct(D2, D2.detected);
        
        [D3, D3u, D3a, D2sub, t_bin, bin_count, syn_wf]=calc_h(D2, HW(kk), params);
        out_file=sprintf('%s/tx_test_case_R=%2.1f_HW=%2.1f_seg%d.h5', out_dir, R(kk), HW(kk), count);
        for kf=1:length(ph_fields);
            this_field=ph_fields{kf};
            h5create(out_file,['/photon/',this_field], size(D2sub.(this_field)),'datatype','double');
            h5write(out_file,['/photon/', this_field], double(D2sub.(this_field)));
        end
        h5create(out_file,'/hist/t_bin', size(t_bin{1}),'datatype','double'); 
        h5write(out_file,'/hist/t_bin', double(t_bin{1}));
        h5create(out_file,'/hist/bin_count', size(bin_count{1}),'datatype','double');
        h5write(out_file,'/hist/bin_count', double(bin_count{1}));
        h5create(out_file,'/h_trunc_corr/median', 1, 'datatype','double');
        h5write(out_file,'/h_trunc_corr/median', D3.med);
        h5create(out_file,'/h_trunc_corr/mean', 1, 'datatype','double');
        h5write(out_file,'/h_trunc_corr/mean', D3.centroid);
        h5create(out_file,'/h_trunc_corr/N_signal_est', 1, 'datatype','double');
        h5write(out_file,'/h_trunc_corr/N_signal_est', D3.N_signal_est);
        h5create(out_file,'/h_trunc_corr/N_noise', 1, 'datatype','double');
        h5write(out_file,'/h_trunc_corr/N_noise', D3.N_noise);
        
        h5create(out_file,'/h_trunc/median', 1, 'datatype','double');
        h5write(out_file,'/h_trunc/median', D3u.med);
        h5create(out_file,'/h_trunc/mean', 1, 'datatype','double');
        h5write(out_file,'/h_trunc/mean', D3u.centroid);
        
        % write out the synthetic WF
        h5create(out_file,'/syn_wf/t', size(syn_wf.t), 'datatype','double');
        h5write(out_file,'/syn_wf/t',  syn_wf.t);
        h5create(out_file,'/syn_wf/P', size(syn_wf.P), 'datatype','double');
        h5write(out_file,'/syn_wf/P',  syn_wf.P);
        h5create(out_file,'/syn_wf/mask', size(syn_wf.mask), 'datatype','double');
        h5write(out_file,'/syn_wf/mask',  double(syn_wf.mask));
        
        
        h5create(out_file,'/params/R', 1, 'datatype','double');
        h5write(out_file,'/params/R', R(kk));
        h5create(out_file,'/params/HW', 1, 'datatype','double');
        h5write(out_file,'/params/HW', HW(kk));
        h5create(out_file,'/params/Hctr', 1,'datatype', 'double')
        h5write(out_file,'/params/Hctr', D3.Hctr);      
    end
end

%-------------------------------------------------------
function S=disp_test_case(dirname)

h5_files=dir([dirname,'/*.h5']);
for k=1:length(h5_files);
    temp=regexp(h5_files(k).name,'tx_test_case_R=(\d+.\d+)_HW=(\d+.\d+)_seg(\d+).h5','tokens');
    R(k)=str2num(temp{1}{1}); 
    HW(k)=str2num(temp{1}{2}); 
    seg(k)=str2num(temp{1}{3});
end
[s, ind]=sort(R+HW/1000+seg/1e6);
h5_files=h5_files(ind);
fields={'params/R', 'params/HW', 'params/Hctr', 'h_trunc/mean','h_trunc/median', 'h_trunc_corr/mean', 'h_trunc_corr/median'};
fprintf(1, '%s\t', 'filename','R','HW','H_ctr','uncorr mean' ,'uncorr median', 'corr_mean','corr median')
fprintf(1,'\n');
for k=1:length(h5_files);
    for kf=1:length(fields);
        S(k, kf)=h5read([dirname,'/', h5_files(k).name],['/',fields{kf}]); 
    end
    fprintf(1, ['%s', repmat('\t%3.3f', [1, kf]),'\n'],  h5_files(k).name, S(k,:));
end
    
uR=unique(S(:,1));
for k=1:length(uR); 
    these=S(:,1)==uR(k); 
    fprintf(1, 'R=%f\tHmean_u=%f\tHmean_c=%f\tHmed_u=%f\tHmed_c=%f\n', uR(k),  mean(S(these, [4 6 5 7]))); 
end


%-------------------------------------------------------
function [D2, params]=make_data(N_pulses, N_chan, roughness, N_per_pulse, WF, BGR)

DEM=0;
x0=zeros(N_pulses,1);
ATM_xmit=ones(size(x0));
params=struct('roughness', roughness,   'sigma_x', 7.5, 'NoiseRate', BGR, 'H_window', 2*max(abs(WF.t))*1.5e8, 'WF', WF, 'c', 3e8, 'N_channels', N_chan, 'N_per_pulse', N_per_pulse,'t_dead', 3.2e-9); 
D2=det_sim(DEM, x0, params, ATM_xmit);
D2.pulse=D2.pulse_num;
D2=rmfield(D2, 'pulse_num');

%-------------------------------------------------------
function [D3, D3u, D3a, D2sub, t_bin, bin_count, syn_wf]=calc_h(D2, HW, params)

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
[D3a.med, D3a.centroid, D3a.count]=deal(NaN(length(segs)-1,1));
D3=D3a;
%[D3.sigma_med, D3.sigma_centroid]=deal(NaN(length(segs)-1,1));
%D3u=D3;
CONSTANTS.dt_hist=25e-12;
if length(segs)==1
    segs(2)=segs(1)+58;
end
[t_bin, bin_count]=deal(cell(length(segs)-1,1));
for k=1:length(segs)-1
    in_bin=D2.pulse>=segs(k) & D2.pulse <=segs(k+1);
    D3a.med(k)=median(D2.h(in_bin));
    D3a.centroid(k)=mean(D2.h(in_bin));
    D3a.count(k)=sum(in_bin);
    
    % take an initial subset of pulses in the segment
    D2sub0=index_struct(D2, in_bin);
    
    % iterate 5x to find the truncated median:
    Hmed=median(D2sub0.h);
    in_HW=abs(D2sub0.h-Hmed)<HW/2;
    for k_it=1:5
        Hmed=median(D2sub0.h(in_HW));
        in_HW=abs(D2sub0.h-Hmed)<HW/2;
    end
    % additional iteration to recenter on the mean
    in_HW=abs(D2sub0.h-mean(D2sub0.h(in_HW)))<HW/2;
    Hmed=median(D2sub0.h(in_HW));
    
    D2sub=index_struct(D2sub0, in_HW);
    
    % calc the uncorrected statistics;
    D3u.med(k)=median(D2sub.h);
    D3u.centroid(k)=mean(D2sub.h);
    [D3.count(k), D3u.count(k)]=deal(length(D2sub.h));
    
    sigma_rx=robust_peak_width_CDF(D2sub.h/(-1.5e8), 58*HW/1.5e8*params.NoiseRate, HW*[-0.5 0.5]/(1.5e8));
    
    % correct the data
    N_bins=ceil((HW/(1.5e8)/CONSTANTS.dt_hist/2))*2;
    t_bin{k}=mean(D2sub.h/(-1.5e8))+(-N_bins/2:N_bins/2)*CONSTANTS.dt_hist;
    %[t_bin{k}, bin_count{k}]=quick_hist(D2sub.h/(-1.5e8), CONSTANTS.dt_hist);
    bin_count{k}=my_histc(D2sub.h/(-1.5e8), t_bin{k});
    %[dM, dCtr]=correct_for_TX_shape(t, z, t_TX, TX, HW, Nrate, N_sig)
    D3.N_noise(k)=58*HW/1.5e8*params.NoiseRate;
    D3.N_signal_est(k)=sum(bin_count{k})-D3.N_noise(k);
    D3.SNR=max(0,D3.N_signal_est(k)./max(1,D3.N_noise(k)));
    [dM, dCtr, syn_wf]=correct_for_TX_shape(params.WF.t, params.WF.p, HW/1.5e8, sigma_rx, D3.SNR);
    %!!!!!!WRONG!!!!!  -- this med correction should be the difference
    %between the TXP WF med and this broadened median.
    D3.med(k)=D3u.med(k)+dM;
    D3.centroid(k)=D3u.centroid(k)+dCtr;
    D3.Hctr=Hmed;
end


%-----------------------------------------------------------
function [bin_ctr, count]=quick_hist(x, dx, x0)

% fast algorithm for generating histograms

if ~exist('x0','var');
    x0=0;
end
xb=sort(round((x-x0)/dx));
bin_num=xb-xb(1)+1;

[uB, ind1]=unique(bin_num,'first');
[uB, ind2]=unique(bin_num,'last');
bin_ctr=(xb(1):xb(end))*dx+x0;
count=zeros(size(bin_ctr));
count(uB)=ind2-ind1+1;




%% NOTE: this code reads the correction data back in and checks it against the correction by hand
if false
    
    %[dM, dCtr]=correct_for_TX_shape(t, z, t_TX, TX, HW, Nrate, N_sig)
    [~, out]=unix(['ls ', out_dir,'/*.h5']); out=strsplit(deblank(out));
    
    t_wf=h5read(out{1},'/hist/t_bin');
    N_wf=h5read(out{1},'/hist/bin_count');
    N_noise=h5read(out{1},'/h_trunc_corr/N_noise');
    N_sig=h5read(out{1},'/h_trunc_corr/N_signal_est');
    HW=h5read(out{1},'/params/HW');
    Nrate=1e7*58;
    [dM, dCtr]=correct_for_TX_shape(t_wf, N_wf, WF.t, WF.p, HW/1.5e8, Nrate, N_sig);
    dM_file=h5read(out{1},'/h_trunc_corr/median')-h5read(out{1},'/h_trunc/median')
    dCtr_file=h5read(out{1},'/h_trunc_corr/mean')-h5read(out{1},'/h_trunc/mean')
    
    t_syn=h5read(out{1},'/syn_wf/t');
    p_syn=h5read(out{1},'/syn_wf/P');
end
