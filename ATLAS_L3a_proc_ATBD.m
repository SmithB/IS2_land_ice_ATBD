function  varargout=ATLAS_L3a_proc_ATBD(varargin)

% wrapper function for ATL03->ATL06 processing.   
if isstruct(varargin{1})
    [varargout{1:nargout}]=feval('ATL03_to_ATL06', varargin{:});
else
    if nargout>0
        [varargout{1:nargout}]=feval( varargin{:});
    else
        feval(varargin{:});
    end
end

%-------------------------------------------------------

function [D3a, dh_hist, LOG]=ATL03_to_ATL06(D2, params, which_seg)

% ATLAS_L3a_proc_ATBD(D2, params, which_L)
% Inputs:
% D2: a 2-element structure array containing two beams worth of L2b
%    ICESat-2 data
%    must contain, at minimum, fields:
%        x_RGT: along-track coordinate
%        h: surface height
%        detected: whether each photon was detected (set to all true if no
%        detector model used
%         ...probably others...
% params: a 2-element (one per beam structure containing metadata
%    about the D2 data and the L2->L3 processing
%      must contain, at minimum, fields:
%        WF: A structure array describing the transmit WF, with fields:
%            t: photon-bin time
%            p: mean power in the bin
%     if params contains a varaible 'SAVELOG' and it is true, the
%       'LOG' output will be generated, and will contain flag values from 
%        different stages of the fitting process.  See below for the 
%        variable 'SAVELOG'
% which_seg : An optional parameter specifying the segment to be processed.
%    If specified, only the segment centered on which_L will be processed.
%    This is useful for plotting and debugging, can be ignored for batch
%    processing.


if isfield(params(1),'SAVELOG') && params(1).SAVELOG
    SAVELOG=true;
else
    SAVELOG=false;
end

CONSTANTS.dt_hist=100e-12;
CONSTANTS.c=299792458; %speed of light (m/s)

if ~isfield(params,'pulses_per_seg'); for kB=1:2; params(kB).pulses_per_seg=57; end;end

if ~isfield(params,'sigma_pulse')
    for k=1:length(params)
        params(k).sigma_pulse=diff(wf_percentile(params(k).WF.t, params(k).WF.p, [0.16 .84]))/2;
    end
end
if isfield(D2,'detected')
    D2_including_undetected=D2;
    for kB=1:2
        D2(kB)=index_struct(D2(kB),D2(kB).detected==1);
    end
end


if isfield(D2,'elev')
    for k=1:length(D2)
        D2(k).h=D2(k).elev;
    end
end

% NEW 10/2015: Throw away PE around gaps in the DEM
for kB=1:2
    if isfield(D2(kB),'z0')
        D2(kB)=index_struct(D2(kB), isfinite(D2(kB).h) & isfinite(D2(kB).z0));
    else
        D2(kB)=index_struct(D2(kB), isfinite(D2(kB).h));
    end
end

% NEW 10/2015: exit if x_RGT is of zero length
x_AT_range=range(real(cat(1, D2.x_RGT)));
if isempty(x_AT_range)
    D3=[]; D3a=[]; dh_hist=[];
    return
end

% new 11/2016: use the segment number if it exists (OTW add)
if ~isfield(D2,'seg_num')
    for kB=1:2
        D2(kB).seg_num=1+floor(D2(kB).x_RGT/20);
    end
end

seg_list=unique(cat(1, D2.seg_num));

if exist('which_seg','var')
    seg_list=seg_list(find(abs(seg_list-which_seg)==min(abs(seg_list-which_seg)), 1, 'first'));
end

% disp('READING L0 from file');
% L0=load('bin_centers.txt');

% NEW: 4/2016:  Add the histogram
dh_hist_full(1:2)=deal(struct('dh', [-4.99:0.01:5]', 'count', uint8(zeros(length(-4.99:.01:5), length(seg_list)))));

% added x_PS_ctr, y_PS_ctr, x_RGT, y_RGT, SNR, exit_iteration
D3_fields={ 'segment_id', 'signal_selection_source', ...
    'signal_selection_status_confident','signal_selection_status_all', 'signal_selection_status_backup', ...
    'h_expected_rms', 'sigma_geo_AT', 'sigma_geo_XT', ...
    'RGT','GT','PT','cycle','orbit_number','seg_count', ...
    'track','beam', 'BGR','h_initial','w_surface_window_initial','w_surface_window_final',  ...
    'dh_fit_dx', 'dh_fit_dy', 'h_robust_spread', 'h_rms','h_mean','sigma_h_mean','h_med', ...
    'sigma_dh_fit_dx', 'fpb_error' ,'fpb_med_corr','fpb_mean_corr', 'med_r_fit', ...
    'N_initial', 'N_noise', 'n_fit_photons','fpb_N_corr', 'sigma_photon_est', ...
    'TX_med_corr','TX_mean_corr', 'h_LI', 'h_LI_sigma', 'z0','m0', 'x_RGT','y_RGT', ...
    'lat_ctr','lon_ctr', ...
    'fpb_med_corr_sigma','fpb_mean_corr_sigma', ...
    'first_seg_pulse','N_seg_pulses','time','SNR','exit_iteration', ...
    'surface_detection_failed','surface_height_invalid','cross_track_slope_invalid' ,'fpb_correction_failed'};

for kf=1:length(D3_fields)
    D3_empty.(D3_fields{kf})=NaN;
end

% NEW 10/2016: initialize ATL06 status and signal_selection_flags
for field={'signal_selection_source', 'signal_selection_status_confident','signal_selection_status_all', 'signal_selection_status_backup' ...
        'surface_detection_failed','surface_height_invalid','cross_track_slope_invalid' ,'fpb_correction_failed'  }
    D3_empty.(field{1})=uint8(0);
end

D3=repmat(D3_empty, [length(seg_list),2]);

% new 4/2016: setup the output version of the dh_hist.
dh_hist.seg_0=find(mod(seg_list, 10)==0);
dh_hist.seg_1=min(dh_hist.seg_0+9, length(seg_list));
dh_hist.dh=dh_hist_full.dh;
dh_hist.count=uint16(zeros(numel(dh_hist.dh), numel(dh_hist.seg_0)));
dh_hist.geo_bin_list=NaN(10, numel(dh_hist.seg_0));
[dh_hist.N_pulses, dh_hist.x_RGT_mean, dh_hist.x_PS_mean, dh_hist.y_PS_mean]=deal(NaN(1, numel(dh_hist.seg_0)));
dh_hist=repmat(dh_hist, 1, 2);


tic
for k0=1:length(seg_list)
    % initial processing: find the intial vertical bin centers
    ybar=NaN(1,2);
    for kB=1:2
        D3(k0, kB).segment_id=seg_list(k0);
        [D3(k0, kB).x_RGT, x0]=deal(seg_list(k0)*20-20);
        AT_els=D2(kB).seg_num==seg_list(k0)-1 | D2(kB).seg_num==seg_list(k0);
        
        if isfield(params(kB),'skip') &&  params(kB).skip
            AT_els=[];
        end
        
        if SAVELOG; LOG(k0, kB).AT_els=AT_els; end
        
        if ~any(AT_els)
            D3(k0, kB).signal_selection_source=3;
            D3(k0, kB).signal_selection_status_confident=3;
            D3(k0, kB).signal_selection_status_all=3;
            D3(k0, kB).signal_selection_status_backup=3;
            continue;
        end
        D2sub_unfilt(kB)=index_struct(D2(kB), AT_els);
        [initial_fit_els, D3(k0, kB).signal_selection_source, D3(k0, kB).signal_selection_status_confident,  D3(k0, kB).signal_selection_status_all]=...
            choose_ground_strategy(D2sub_unfilt(kB));
 
        if D3(k0, kB).signal_selection_source < 1  % enough confident PE found to determine slope and height of the bin
            initial_Hwin_min=3;
        else %
            initial_Hwin_min=10;
        end
        
        if D3(k0, kB).signal_selection_source <=1  % if we are using ATL03 flagged PE
            [initial_fit_els, D3(k0, kB).w_surface_window_initial, D3(k0, kB).h_initial, ~]=...
                initial_at_fit(D2sub_unfilt(kB).x_RGT, D2sub_unfilt(kB).h, initial_fit_els, x0, median(D2sub_unfilt(kB).BGR), initial_Hwin_min, params(kB));
            if D3(k0, kB).signal_selection_source==1
                D2sub(kB)=index_struct(D2sub_unfilt(kB), initial_fit_els);
                initial_fit_els=true(sum(initial_fit_els),1);
            else
                D2sub(kB)=D2sub_unfilt(kB);
            end
        end
        
        if D3(k0, kB).signal_selection_source ==2 &&  D3(k0, kB).signal_selection_status_backup ~=3  % not enough confident or padded PE have been found, fall back to alternate strategies
            [initial_fit_els, D3(k0, kB).signal_selection_source, D3(k0, kB).signal_selection_status_backup]=...
                backup_signal_finding_strategy(D2sub_unfilt(kB), D2(kB), seg_list(k0), 20);
            D2sub(kB)=index_struct(D2sub_unfilt(kB), initial_fit_els);
            initial_fit_els=true(size(D2sub(kB).h));
            D3(k0, kB).w_surface_window_initial=10;
            D3(k0, kB).h_initial=mean(D2sub(kB).h(initial_fit_els));
        end
        if SAVELOG; LOG(k0, kB).initial_fit_els=initial_fit_els; end
                
        if D3(k0, kB).signal_selection_source<3 % go ahead if we have initial PE
            if SAVELOG
                LS_fit_options=struct( 'Nsigma', 3, 'Hwin_min', 3,'SAVELOG', true);
            else
                LS_fit_options=struct( 'Nsigma', 3, 'Hwin_min', 3);
            end
            max_ground_bin_iterations=100;
            if SAVELOG
                [D3(k0, kB), r, els, LOG(k0, kB).LS_fit]=ATLAS_LS_fit(D2sub(kB), D3(k0, kB).x_RGT, initial_fit_els, max_ground_bin_iterations, params(kB), D3(k0, kB), LS_fit_options);
            else
                [D3(k0, kB), r, els]=ATLAS_LS_fit(D2sub(kB), D3(k0, kB).x_RGT, initial_fit_els, max_ground_bin_iterations, params(kB), D3(k0, kB), LS_fit_options);
            end
             
            if SAVELOG
                LOG(k0, kB).els_after_LS_fit=els;
                LOG(k0, kB).r_after_LS_fit=r;
            end
            
            D3(k0, kB).med_r_fit=median(r);
            D3(k0, kB).n_fit_photons=length(r);
            
            % make the time histogram:
            clear temp;
            temp.time=-2/CONSTANTS.c*r;
            t_WF=((min(temp.time)-2.5*CONSTANTS.dt_hist):CONSTANTS.dt_hist:(max(temp.time)+params(kB).t_dead))';
            N_WF=my_histc(  temp.time, t_WF);
            
            if ~isfield(params,'N_channels'); [params(:).N_channels]=deal(params(:).N_det); end;
            % Changed 3/24/2016
            [D3(k0, kB).fpb_med_corr, D3(k0, kB).fpb_mean_corr, D3(k0, kB).fpb_N_corr, N_WF_corr, D3(k0, kB).fpb_med_corr_sigma, D3(k0, kB).fpb_mean_corr_sigma, minGain]=...
                fpb_corr(t_WF, N_WF, params(kB).N_channels, 57, params(kB).t_dead,  CONSTANTS.dt_hist);
            % check if gain correction is valid
            if minGain < 1/(2*params(kB).N_channels)
                D3(k0, kB).fpb_correction_failed=true;
            end
            if D3(k0, kB).signal_selection_source>=2
                D3(k0, kB).surface_detection_failed=1;
            end
            
            D3(k0, kB).N_noise=median(D2sub(kB).BGR)*D3(k0, kB).w_surface_window_final/(CONSTANTS.c/2)*57;
            %sigma_hat_robust=robust_peak_width_from_hist(t_WF, N_WF_corr, D3(k0, kB).N_noise, D3(k0, kB).w_surface_window_final*[-0.5 0.5]/(CONSTANTS.c/2));
            %[D3(k0, kB).TX_med_corr, D3(k0, kB).TX_mean_corr]=correct_for_TX_shape(sigma_hat_robust,[],  params(kB).WF.t, params(kB).WF.p, D3(k0, kB).w_surface_window_final/(1.5e8));
            if isfield(params(kB),'WF')
                [D3(k0, kB).TX_med_corr, D3(k0, kB).TX_mean_corr]=correct_for_TX_shape(t_WF, N_WF_corr, params(kB).WF.t, params(kB).WF.p, D3(k0, kB).w_surface_window_final/(1.5e8),...
                    median(D2sub(kB).BGR)*57, max(10, sum(N_WF_corr)-D3(k0, kB).N_noise));
            else
                D3(k0, kB).TX_med_corr=0;
            end
            D3(k0, kB).h_LI=D3(k0, kB).h_mean + D3(k0, kB).fpb_med_corr + D3(k0, kB).TX_med_corr;
            
            if false
                % extra code to estimate the true fpb error
                GG=[ones(size(D2sub_all(kB).h)), real(D2sub_all(kB).x_RGT)-D3(k0, kB).x_RGT];
                els_all=abs(D2sub_all(kB).h-GG*[D3(k0, kB).h_mean; D3(k0, kB).dh_fit_dx]) < D3(k0, kB).w_surface_window_final/2;
                m_all=GG(els_all,:)\D2sub_all(kB).h(els_all);
                %fprintf(1, 'widths are: [sigma_hat_robust, sigma_hat_robust_all]=\t[%4.2d %4.2d], slope spreading is %3.2d\n', sigma_hat_robust*1.5e8, sigma_hat_robust_all, std(D2sub(kB).zground(D2sub(kB).SigNoise==1)-D2sub(kB).z0(D2sub(kB).SigNoise==1)));
                D3(k0, kB).fpb_error=D3(k0, kB).h_mean-m_all(1);
            end
            % now report the across-track coordinates for the segment
            % NEW 10/2015:  Changed order of operations (x_RGT is calculated
            % earlier)
            if isfield(D2sub,'y_RGT')
                D3(k0, kB).y_RGT=median(D2sub(kB).y_RGT);
            else
                D3(k0, kB).y_RGT=median(imag(D2sub(kB).x_RGT));
            end
            temp=unique(D2sub(kB).beam(isfinite(D2sub(kB).beam)));
            D3(k0, kB).beam=temp(1);
            temp=unique(D2sub(kB).track(isfinite(D2sub(kB).track)));
            D3(k0, kB).track=temp(1);
            D3(k0, kB).first_seg_pulse=min(D2sub(kB).pulse_num);
            D3(k0, kB).N_seg_pulses=diff(range(D2sub(kB).pulse_num))+1;
            D3(k0, kB).time=median(D2sub(kB).time);
            
            % NEW 10/2015:  Use regression to calculate reference center
            % parameters
            % regress WRT along-track distance to get parameters for segment
            % center
            %'x_PS_ctr', 'y_PS_ctr', 'lat_ctr','lon_ctr'
            Ginv_AT=[ones(size(D2sub(kB).x_RGT(els))) D2sub(kB).x_RGT(els)-D3(k0, kB).x_RGT]\eye(sum(els));
            Ginv_AT=Ginv_AT(1,:);
            %             D3(k0, kB).x_PS_ctr=Ginv_AT*D2sub(kB).x0(els);
            %             D3(k0, kB).y_PS_ctr=Ginv_AT*D2sub(kB).y0(els);
            D3(k0, kB).lat_ctr=Ginv_AT*D2sub(kB).lat(els);
            D3(k0, kB).lon_ctr=Ginv_AT*D2sub(kB).lon(els);
            
            
            if isnan(D3(k0, kB).h_LI)
                D3(k0, kB).surface_height_invalid=true;
            end
            
            ybar(kB)=median(D2sub(kB).y_RGT);
        end
    end
    % calculate dh_fit_dy:
    dh_fit_dy=(D3(k0,2).h_LI-D3(k0, 1).h_LI)/(diff(ybar));
    for kB=1:2
        D3(k0, kB).dh_fit_dy=dh_fit_dy;
    end
    % error propagation steps, calculated after final dh_dy is set
    
    if ~isfinite(dh_fit_dy); dh_fit_dy=0; end
    for kB=1:2
        if ~isfinite(D3(k0, kB).h_mean); continue; end
        N_sig=max(0,D3(k0, kB).n_fit_photons-D3(k0, kB).N_noise);
        sigma_signal=sqrt((dh_fit_dy.^2+D3(k0, kB).dh_fit_dx.^2)*params(kB).sigma_x.^2 + (params(kB).sigma_pulse*1.5e8).^2);
        D3(k0, kB).sigma_photon_est=sqrt((D3(k0, kB).N_noise*(D3(k0, kB).w_surface_window_final*.287).^2+N_sig*sigma_signal.^2)/D3(k0, kB).n_fit_photons);
        sigma_per_photon=max(D3(k0, kB).sigma_photon_est, D3(k0, kB).h_robust_spread);
        D3(k0, kB).sigma_h_mean=D3(k0, kB).sigma_h_mean*sigma_per_photon;
        D3(k0, kB).sigma_dh_fit_dx=D3(k0, kB).sigma_dh_fit_dx*sigma_per_photon;
        D3(k0, kB).h_LI_sigma=max(D3(k0, kB).sigma_h_mean, D3(k0, kB).fpb_med_corr_sigma);
    end
    
    for kB=1:2
        if isfield(params,'RGT')
            D3(k0, kB).RGT=params(kB).RGT;
            D3(k0, kB).GT=params(kB).GT;
            D3(k0, kB).PT=params(kB).PT;
            D3(k0, kB).cycle=params(kB).cycle;
            D3(k0, kB).orbit_number=params(kB).orbit_number;
        end
    end
    % new, 4/2016: calculate residual histogram
    for kB=1:2
        if ~isfinite(D3(k0, kB).h_mean); continue; end
        xx=D2sub_unfilt(kB).x_RGT-D3(k0, kB).x_RGT;
        zz=D2sub_unfilt(kB).h;
        zz=zz(abs(xx)<10);  % select elements from the non-overlapping parts of the segment
        xx=xx(abs(xx)<10);
        if isempty(xx); continue; end
        G=[ones(size(xx)), xx];
        r=zz-G*[D3(k0, kB).h_mean; D3(k0, kB).dh_fit_dx];
        dh_hist_full(kB).count(:, k0)=uint8(my_histc(r, dh_hist_full(kB).dh));
    end
    
    if mod(k0, 250)==0
        disp([num2str(k0) ' out of ' num2str(length(seg_list))]);
        clear D3a;
        f=fieldnames(D3); for kF=1:length(f); for kB=1:2; D3a.(f{kF})(:,kB)=cat(1, D3(:, kB).(f{kF})); end; end
        r=[D3a.h_LI-D3a.z0];
        fprintf(1, 'mean r =%3.2d, std r =%3.2d\n', [mean(r(isfinite(r))), std(r(isfinite(r)))]);
        %         if exist('save_file','var');
        %             save(save_file, 'D3');
        %         end
    end
end

if isempty(D3)
    D3a=[];
    return
end

%  new 4/2016 : collect the D3 to D3a:
clear D3a;
ff=fieldnames(D3);
for kF=1:length(ff)
    for kB=1:2
        D3a.(ff{kF})(:,kB)=cat(1, D3(:, kB).(ff{kF}));
    end;
end

% new 4/2016: stack the dh_histogram by 10.
for kB=1:2
    for ks=1:length(dh_hist(kB).seg_0)
        els=dh_hist(kB).seg_0(ks):dh_hist(kB).seg_1(ks);
        els=els(D3a.surface_height_invalid(els, kB)==0   & D3a.SNR(els, kB) > 0.5);
        if ~isempty(els)
            dh_hist(kB).segment_id_list(1:length(els), ks)=D3a.seg_count(els, kB);
            dh_hist(kB).count(:, ks)=uint16(sum(dh_hist_full(kB).count(:, els), 2));
            dh_hist(kB).x_RGT_mean(ks)=mean(D3a.x_RGT(els, kB));
            dh_hist(kB).lat_mean(ks)=mean(D3a.lat_ctr(els, kB));
            dh_hist(kB).lon_mean(ks)=mean(D3a.lon_ctr(els, kB));
            dh_hist(kB).N_pulses(ks)=sum(D3a.N_seg_pulses(els, kB))/2;
        end
    end
end

%--------------------------------------------------------------------------
function [initial_els, w_surface_window, h_window_ctr, AT_slope]=initial_at_fit(x_AT, h, initial_els, x0, BGR, W_min, params)

G=[ ones(size(x_AT)), x_AT-x0];
m=G(initial_els,:)\h(initial_els);
r=h(:)-G*m;
Noise_Ph_per_m=params.pulses_per_seg*BGR/(params.c/2);
H_win=diff(range(r(initial_els)));
sigma_r=robust_peak_width_CDF(r(initial_els), Noise_Ph_per_m*H_win, [0 H_win]-H_win/2);
sigma_expected=sqrt((params.sigma_x*m(2)).^2+ (params.sigma_pulse*params.c/2).^2);
w_surface_window=max(W_min, 6*max(sigma_r, sigma_expected));
h_window_ctr=m(1); AT_slope=m(2);
% pad the surface window if needed
initial_els=initial_els | abs(r) < w_surface_window/2;

%---------------------------------------------------------------------------------------
function [D3, r0, els, LOG]=ATLAS_LS_fit(D2, L0, initial_els,  N_it, params, D3, options)

if isfield(options,'SAVELOG')
    SAVELOG=true;
else
    SAVELOG=false;
end
if ~exist('options','var')
    options=struct('Nsigma', 3, 'Hwin_min', 3);
end

c2=3e8/2;
if isfield(D2,'elev')
    D2.h=D2.elev;
end

m=[]; sub1=[]; r0=[];

if SAVELOG
    LOG.iterations=struct('els', [], 'sigma_r', [], 'sigma_expected', [], 'W_win', [], 'h_ctr', [], 'dhdx', []);
end

% fit an along-track polynomial
%s_ctr=mean(range(real(D2.x_LC)));
ds=real(D2.x_RGT)-L0;
els=initial_els;

G=[ones(size(ds(:))), ds(:)];
m=G(els,:)\D2.h(els);

D3.N_initial=sum(els);
if SAVELOG; LOG.G=G; end
% iterate to reduce residuals
Noise_Ph_per_m=diff(range(D2.pulse_num))*median(D2.BGR)/c2;
H_win=diff(range(D2.h(els)));
for k=1:N_it
    m_last=m;
    m=G(els,:)\D2.h(els);
    
    r_all=D2.h-G*m;
    r0=r_all(els);
    sigma_r=robust_peak_width_CDF(r0, Noise_Ph_per_m*H_win, [0 H_win]-H_win/2);
    sigma_expected=sqrt((c2*params.sigma_pulse).^2+params.sigma_x.^2*(m(2).^2));
    
    els_last=els;
    SNR_last=sum(els)/(H_win*Noise_Ph_per_m);
    H_win_last=H_win;
    H_win=max([2*[sigma_expected, sigma_r]*options.Nsigma, 0.75*H_win_last, options.Hwin_min]);
    els=abs(r_all ) < H_win/2;
    
    if SAVELOG
        LOG.iterations(k).els=els;
        LOG.iterations(k).sigma_r=sigma_r;
        LOG.iterations(k).sigma_expected=sigma_expected;
        LOG.iterations(k).W_win=H_win;
        LOG.iterations(k).h_ctr=m(1);
        LOG.iterations(k).dhdx=m(2);
    end
    
    SNR=sum(els)/(H_win*Noise_Ph_per_m);
    
    % The last || is experimental-- does it help to quit iterating if SNR
    % decreases?
    % A: No.  Not helpful.
    if sum(els) < 10 || diff(range(D2.x_RGT(els))) < 20  %|| SNR < (SNR_last-1/(H_win*Noise_Ph_per_m));
        m=m_last;
        H_win=H_win_last;
        els=els_last;
        break
    end
    
    if sum(abs(els_last-els))==0  % no change from last iteration, or we've converged.
        break
    end
    %find(els_last-els)
end

% plotting command:
%clf; plot(real(D2.x0), D2.h,'.'); hold on; plot(real(D2.x0(els)), D2.h(els),'ro'), plot(real(D2.x0(els)), G(els,:)*m,'g')

D3.h_mean=m(1);
D3.dh_fit_dx=m(2);
D3.BGR=median(D2.BGR);
r0=r_all(els);
D3.w_surface_window_final=H_win;

D3.h_robust_spread=iqr(r0)/2;
D3.h_rms=std(r0);
D3.h_med=m(1)+median(r0);

% NEW 10/2015: Calculate the noise estimate and the SNR
Noise_est=H_win*Noise_Ph_per_m-length(unique(D2.pulse_num))*median(D2.BGR)/c2;
D3.SNR=(sum(els)-Noise_est)./Noise_est;

% NEW 10/2015: Report the exit iteration
D3.exit_iteration=k;
if N_it==1
    D3.exit_iteration=0;
end

G1=G(els,:);
Ginv=(G1'*G1)\G1';
CovMat=Ginv*Ginv';

if SAVELOG
    LOG.CovMat=CovMat;
end

% assign unscaled error estimate
D3.sigma_h_mean=sqrt(CovMat(1,1));
D3.sigma_dh_fit_dx=sqrt(CovMat(2,2));
if isfield(D2,'z0')
    % true height: sample surface every meter, fit with slope
    ds0=round_to(ds, 1); [~, ds0_ind]=unique(ds0);
    G=[ones(length(ds0_ind),1), ds(ds0_ind)];
    m=G\D2.z0(ds0_ind);
    D3.z0=m(1); D3.m0=m(2);
else
    D3.z0=NaN; D3.m0=NaN;
end

% NEW 10/2015: choose_ground_strategy scheme
%---------------------
function [els, signal_selection_source, sss_confident, sss_all]=choose_ground_strategy(D2)

[signal_selection_source, sss_confident, sss_all]=deal(0);

% check if there are enough flagged (not padded) photons
els=D2.ph_class > 1 &  isfinite(D2.h) ;

if ~any(els) || diff(range(D2.x_RGT(els))) < 20
    % not enough space between first and last PE
    sss_confident=1;
    signal_selection_source=1;
end
if  ~any(els) || sum(els) <10
    % not enough flagged PE
    if sss_confident==1
        sss_confident=3;
    else
        sss_confident=2;
    end
    signal_selection_source=1;
end
% if all good, return
if signal_selection_source==0
    return
end

% now check if there are enough PE including the pad
els=D2.ph_class >= 1 &  isfinite(D2.h);
if  ~any(els) || diff(range(D2.x_RGT(els))) < 20
    % along-track spread of PE too small
    signal_selection_source=2;
    sss_all=1;
end
if sum(els) <10
    % number of selected PE too small
    if sss_all==1
        sss_all=3;
    else
        sss_all=2;
    end
    signal_selection_source=2;
end

%-------------------------------------------
function [selected_PE, signal_selection_source, signal_selection_status]=backup_signal_finding_strategy(D2_local, D2_all, seg_num, W)

found=D2_local.ph_class > 1 &   isfinite(D2_local.h);
if any(found)
    % if we have any selected PE, center the window on them, see if this
    % gives us a valid window (of any quality
    Hmean=mean(D2_local.h(found));
    els=abs(D2_local.h-Hmean)<5;
    if sum(els) > 10 && diff(range(D2_local.x_RGT(els)))>20
        selected_PE=els;
        signal_selection_source=2;
        signal_selection_status=0;
        return
    end
end

% nothing found, or the PE around the found PE were not usable.

these=D2_all.seg_num >=  seg_num-2 & D2_all.seg_num <= seg_num+1;
if sum(these)<10
    selected_PE=[];
    signal_selection_source=3;
    signal_selection_status=4;
    return
end

% use the histogram strategy
h=D2_all.h(these);
bins=(floor(min(h))+0.25):0.5:ceil(max(h));
count=my_histc(h, bins);
% lazy programming: the count for bin +-W is found by convolving with a
% boxcar 2W high
C1=conv(count(:), ones(2*W, 1),'same');

% procede if more than 15 PE are in the best bin
if max(C1) > 20   
    % select the bins that are not significantly differnt from the peak
    z0r=range(bins(C1>max(C1)-sqrt(max(C1))));
    selected_PE=D2_local.h >= z0r(1)-W/2 & D2_local.h <= z0r(2)+W/2;
else
    selected_PE =[];
end
if sum(selected_PE) > 10 &&  diff(range(D2_local.x_RGT(selected_PE)))>20
    signal_selection_source=2;
    signal_selection_status=1;
    return
end
% If we haven't returned yet, the signal selection has failed, figure out
% why
signal_selection_source=3;
signal_selection_status=0;
if diff(range(D2_local.x_RGT(selected_PE)))<20
    signal_selection_status=2;
end
if sum(selected_PE) < 10
    if signal_selection_status==2
        signal_selection_status=4;
    else
        signal_selection_status=3;
    end
end





