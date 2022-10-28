function WF=proc_TEP(in_file, out_name, WF)

% read in the waveform estimate, if none is provided
if ~exist('WF','var');
    if exist('in_file','var') && ~isempty(in_file);
        D=xlsread(in_file, 'B2:c1001');
    else
        D=xlsread('TEPforBen-10to33ns.xlsx','B2:C462');
    end
    
    WF.t=D(:,1);
    WF.p=D(:,2);
end

% start with the part of the signal between times 1e-8 and 3.3e-8
% this is the TEP spec
els=WF.t>=10*1e-9 & WF.t < 33e-9;
WF.t=WF.t(els); WF.p=WF.p(els);

% Assume that within this window, the signal is between 5 and 10 ns WRT the
% start
sig_samps=(WF.t-min(WF.t) > 0.5e-8 & max(WF.t)-WF.t > 1.e-8);
noise_samps=~sig_samps;
% initial values
t_ctr=mean(WF.t);
sigma=mean(WF.t)/2;

% iterate ten times to select the range of samples within 3 sigma of the
% pulse median
for k=1:10
    N=mean(WF.p(noise_samps));
    P1=WF.p-N;
    %P1=max(eps*10,WF.p-N);
    t_ctr_last=t_ctr;
    t_ctr=sum(WF.t(sig_samps).*P1(sig_samps))/sum(P1(sig_samps));
    sigma_last=sigma;
    sigma=(wf_percentile(WF.t(sig_samps), P1(sig_samps), 0.84)-wf_percentile(WF.t(sig_samps), P1(sig_samps), 0.16))/2;
    sig_samps=abs(WF.t-t_ctr)<6*sigma;
    if all(abs([t_ctr_last, sigma_last]-[t_ctr, sigma])<1e-13);
        break
    end
end

% center the waveform on its centroid
WF.t=WF.t-t_ctr;
% set bins that are less than (a small number) to (a small number)
WF.p=max(P1, max(P1)*1e-13);
% set bins outside the 3-sigma window to (a small number)
WF.p(~sig_samps)=max(P1)*1e-13;

% write out the waveform estimate, if desired.
if exist('out_name','var') && ~isempty(out_name);
    h5create(out_name, '/time', length(WF.t),'datatype','double');
    h5write(out_name, '/time',double(WF.t));
     
    h5create(out_name,'/power',length(WF.t), 'datatype','double');
    h5write(out_name, '/power', double(WF.p));
end




