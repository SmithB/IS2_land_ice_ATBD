function D2=det_sim(DEM, x0, params, ATM_xmit, CCR)
% the top-level function is a wrapper for the function sub_sim
%
% inputs : DEM: a struct with fields: x: column x coordinate
%                                     y: row y coordinate
%                                     z: surface height
%
%          x0: pulse center coordinates
% params: a sturct with simulation parameter fields:
%      'roughness', surface roughness, meters
%      'sigma_x',  footprint width, meters
%      'NoiseRate', noise rate, in Hz
%      'H_window', size of surface window, meters
%      'WF', struct giving transmit-pulse power as a function of time, with fields t and p
%      'c',  speed of light, m/s
%      'N_channels', number of channels in the detector
%      'N_per_pulse', mean number of detected photons in the return pulse
%      't_dead', 3.2e-9 dead time, in seconds
%
% ATM_xmit: optional value giving the atmospheric attenuation (1 means no attenuation)
%
% CCR: Corner-cube reflector list.  For each, give 
%      'x' - pulse-center location
%      'h' - height
%      'N' - Expected number of PE for a return at the beam center
%
% outputs: D2: struct array giving ATL03-type output data 
%    

if ~exist('ATM_xmit','var')
    ATM_xmit=ones(size(x0));
end

% break the simulation into bite-sized pieces (OTW the inefficiency in appending data becomes the driving cost)
sub_ind=1:1000:length(x0);
if sub_ind(end)<length(x0)
    sub_ind(end+1)=length(x0)+1;
end

if isfield(params,'DEBUG') && params.DEBUG
    % only do the first 10000
    sub_ind=sub_ind(1:min(10, length(sub_ind)));
end

% run the simluation for each piece
for k=1:length(sub_ind)-1;
    x0_sub=x0(sub_ind(k):sub_ind(k+1)-1);
    if isstruct(DEM)
        this_W=200;
        DEM_sub=subset_image(DEM, range(real(x0_sub))+[-1 1]*this_W, range(imag(x0_sub))+[-1 1]*this_W);
        while this_W < 2e4 && (isempty(DEM_sub.x) || isempty(DEM_sub.y))
            this_W=this_W*2;
            DEM_sub=subset_image(DEM, range(real(x0_sub))+[-1 1]*this_W, range(imag(x0_sub))+[-1 1]*this_W);
        end
    else
        DEM_sub=DEM;
    end
    D(k)=sub_sim(DEM_sub, x0(sub_ind(k):sub_ind(k+1)-1), params, ATM_xmit(sub_ind(k):sub_ind(k+1)-1));
    D(k).pulse_num=D(k).pulse_num+sub_ind(k)-1;
end

% catenate the output
f=fieldnames(D);
for k=1:length(f)
    D2.(f{k})=cat(1, D.(f{k}));
end


%-----------------------------------------------
function D=sub_sim(DEM, x0, params, ATM_xmit)
 
D0=struct('xground',[],'pulse_num', [], 'x0', [],'SigNoise', [], 'z0', []);
for k=1:length(x0);
    % Poisson-distributed N
    N_sig=poisson_rv(params.N_per_pulse*ATM_xmit(k), 1);  
    % random footprint location
    D0(k).xground=x0(k)+(randn(N_sig,1)+1i*randn(N_sig,1))*params.sigma_x;    
    D0(k).pulse_num=k+zeros(N_sig,1);
    D0(k).x0=x0(k)+zeros(N_sig,1);
    D0(k).SigNoise=true(N_sig,1);
    
    N_noise=params.NoiseRate*params.H_window/(params.c/2);
    NoiseCount=poisson_rv(N_noise,1);
    if NoiseCount>0;
        D0(k).xground=[D0(k).xground; x0(k)+zeros(NoiseCount,1)];
        D0(k).pulse_num=[D0(k).pulse_num; k+zeros(NoiseCount,1)];
        D0(k).x0=[D0(k).x0; x0(k)+zeros(NoiseCount,1)];
        D0(k).SigNoise=[D0(k).SigNoise; false(NoiseCount,1)];
    end
end

D.x0=cat(1, D0.x0);
D.xground=cat(1, D0.xground);
D.pulse_num=cat(1, D0.pulse_num);
D.SigNoise=cat(1, D0.SigNoise);

if isfield(DEM,'xsub');
    D.zground=interp2(DEM.x, DEM.y, double(DEM.z), real(D.xground), imag(D.xground),'*linear');
    
    NoiseEls=D.SigNoise==0;
    NoiseCtr=interp2(DEM.xsub, DEM.ysub, double(DEM.zsub), real(D.xground(NoiseEls)), imag(D.xground(NoiseEls)),'*linear');
    NoiseGround=D.zground(NoiseEls);
    bad=NoiseCtr-NoiseGround>(params.H_window/2-5);
    NoiseCtr(bad)=NoiseGround(bad)+params.H_window/2-5;
    bad=NoiseCtr-NoiseGround<(-params.H_window/2+5);
    NoiseCtr(bad)=NoiseGround(bad)-params.H_window/2+5;    
    D.zground(NoiseEls)=NoiseCtr+(rand(sum(NoiseEls),1)-0.5)*params.H_window;
else
    if ~isstruct(DEM);
        D.zground=zeros(size(D.xground))+DEM;
    else
        D.zground=interp2(DEM.x, DEM.y, double(DEM.z), real(D.xground), imag(D.xground),'*linear');
    end
    D.zground(D.SigNoise==0)=D.zground(D.SigNoise==0)+(rand(sum(D.SigNoise==0),1)-0.5)*params.H_window;
end
D.zground=D.zground(:);
if ~isstruct(DEM);
    D.z0=zeros(size(D.x0))+DEM;
else
    D.z0=interp2(DEM.x, DEM.y, double(DEM.z), real(D.x0), imag(D.x0),'*linear');
end

if isfield(params, 'WF');
    D.t_ph=-D.zground/(params.c/2)+random_wf_sim(params.WF.t, params.WF.p, length(D.zground));
else
    D.t_ph=-D.zground/(params.c/2)+randn(size(D.zground))*params.sigma_pulse;
end


if isfield(params,'roughness')
    D.t_ph=D.t_ph+randn(size(D.t_ph))*params.roughness*2/params.c;
end

D.h=-D.t_ph*params.c/2;


[~, ind]=sort(D.t_ph + 1e-4*D.pulse_num);
D=index_struct(D, ind);

D.BGR=zeros(size(D.x0))+params.NoiseRate;

[D.detected, D.channel]=detect_with_deadtime(D.t_ph, D.pulse_num, params.N_channels, params.t_dead, length(x0));


% need fields : x_LC, elev, BGR, pulse_num



%----------------------------------------------------
function [hit_all, det, dead_all] = detect_with_deadtime(t, pulse,  N_det, t_dead, N_pulses, t_WF)
% NOTE: t and pulse must be sorted by pulse, then by time.
DS=ceil(rand(size(t))*N_det)+1i*pulse;
u_DS=unique(DS);

if exist('t_WF','var');
    dead_all=zeros(length(t_WF), N_pulses*N_det);
    dt=t_WF(2)-t_WF(1);
    N_dead_bins=ceil(t_dead/dt);
    calc_gain=true;
else
    calc_gain=false;
end

hit_all=false(size(t));
for k=1:length(u_DS);
%     if mod(k, 1e4)==0
%         fprintf(1, '%d out of %d, last_dt = %f\n', k, length(u_DS), toc); 
%         tic
%     end
    these=DS==u_DS(k);
    if ~any(these); 
        continue; 
    end
    % ti are the photons for this detector
    ti=t(these);
    hit=false(size(ti));
    missed=false(size(ti));
    next_first=1;
    hit(next_first)=true;
    
    if calc_gain
        dead=zeros(size(t_WF));
        dead(t_WF>ti(1) & t_WF < ti(1)+t_dead) = true;
    end
    %disp(k)
    while ~(hit(end) || missed(end))
        first=next_first;
        hit(first)=true;
        in_dt=ti>ti(first) & ti < ti(first)+t_dead;
        if any(in_dt)
            % extend the dead time if there are photon hits while the counter is dead
            count=0;
            while sum(in_dt) > count
                count=sum(in_dt);
                in_dt=ti > ti(first) & ti < max(ti(in_dt))+t_dead;
            end
            missed(in_dt) = true;
            last=find(in_dt, 1, 'last');
            if calc_gain
                dead(t_WF>ti(first) & t_WF < ti(last)+t_dead) = true;
            end
            next_first=last+1;
        else
            if calc_gain
                dead(t_WF>ti(first) & t_WF < ti(first)+t_dead) = true;
            end
            next_first=first+1;
        end
    end
    if calc_gain
        dead_all(:, real(u_DS(k))+imag(u_DS(k))*N_det)=dead(:);
    end
    hit_all(these)=hit;
end
det=real(DS);

%_____________________________________________
function D=subset_image(D, XR, YR)

rows=D.y >=YR(1) & D.y <= YR(2);
cols=D.x >= XR(1) & D.x <= XR(2);

D.x=D.x(cols); 
D.y=D.y(rows); 
D.z=D.z(rows, cols);

