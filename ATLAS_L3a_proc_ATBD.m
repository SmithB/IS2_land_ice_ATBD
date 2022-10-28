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
function [D3a, dh_hist, LOG]=ATL03_to_ATL06(D2, params, which_seg, SNR_F_table, dist_for_segment)

% ATLAS_L3a_proc_ATBD(D2, params, which_seg)
% Inputs:
% D2: a 2-element structure array containing two beams worth of L2b
%    ICESat-2 data
%    must contain, at minimum, fields:
%        x_RGT: along-track coordinate WRT the reference ground track
%        h: surface height
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
% SNR_F_table: a look-up table giving the probability that a segment with
%     the observed SNR would be found for a given background rate and
%     initial window size
% dist_for_segment:  a sparse array giving the along-track distance for the
%     segment_id numbers in the input file

if isfield(params(1),'SAVELOG') && params(1).SAVELOG
    SAVELOG=true;
else
    SAVELOG=false;
end

%CONSTANTS.dt_hist=5e-12;
% NOTE: CHANGED 4/24/2019 to match ASAS
CONSTANTS.dt_hist=2.5e-11;

CONSTANTS.c=299792458; %speed of light (m/s)
CONSTANTS.residual_histogram_spacing=0.01;

for kB=1:2
    tt=min(params(kB).WF.t):CONSTANTS.dt_hist:max(params(kB).WF.t);
    params(kB).WF=struct('t', tt, 'p', interp1(params(kB).WF.t, params(kB).WF.p, tt));
end

%if ~isfield(params,'pulses_per_seg'); for kB=1:2; params(kB).pulses_per_seg=57; end;end

if ~isfield(params,'sigma_pulse')
    for k=1:length(params)
        params(k).sigma_pulse=diff(wf_percentile(params(k).WF.t, params(k).WF.p, [0.16 .84]))/2;
    end
end

% in some processing schemes, detected and undetected photons are provided,
% and some are marked as undetected.
if isfield(D2,'detected')
    for kB=1:2
        D2(kB)=index_struct(D2(kB),D2(kB).detected==1);
    end
end

% removed conversion from D2.h_ph and D2.elev to D2.h (must use h_ph)

% removed conversion from ph_class to signal_conf_ph

% NEW 10/2015: Throw away PE around gaps in the DEM
for kB=1:2
    if isfield(D2(kB),'z0')
        D2(kB)=index_struct(D2(kB), isfinite(D2(kB).h_ph) & isfinite(D2(kB).z0));
    else
        D2(kB)=index_struct(D2(kB), isfinite(D2(kB).h_ph));
    end
end

% NEW 10/2015: exit if x_RGT is of zero length
x_AT_range=range(real(cat(1, D2.x_RGT)));
if isempty(x_AT_range)
    D3=[]; D3a=[]; dh_hist=[];
    return
end

%Removed checks on prameters alternative to segment_id. The script should crash if nothing is provided 
    
if ~exist('dist_for_segment','var')
    uSeg=unique(cat(1, D2.segment_id));
    uSeg=uSeg(isfinite(uSeg));
    for kB=1:2
         dist_for_segment{kB}=sparse(uSeg, ones(size(uSeg)), (uSeg-1)*20);
    end
end

seg_list=unique(cat(1, D2.segment_id));
seg_list=seg_list(seg_list < min(length(dist_for_segment{1}), length(dist_for_segment{2})));

if exist('which_seg','var') && ~isempty(which_seg)
    seg_list=unique(seg_list(ismember(seg_list, which_seg)));
end

% NEW: 4/2016:  Add the histogram
dh_hist_full(1:2)=deal(struct('dh', [-4.99:0.01:5]', 'count', uint8(zeros(length(-4.99:CONSTANTS.residual_histogram_spacing:5), length(seg_list)))));

% added x_PS_ctr, y_PS_ctr, x_RGT, y_RGT, SNR, exit_iteration
D3_fields={  'ATL06_quality_summary', 'signal_selection_source', ...
    'signal_selection_status_confident','signal_selection_status_all', 'signal_selection_status_backup', ...
    'h_expected_rms', 'sigma_geo_AT', 'sigma_geo_XT', ...
    'RGT','PT','cycle','seg_count', ...
    'track','beam', 'BGR','h_initial','w_surface_window_initial','w_surface_window_final',  ...
    'dh_fit_dx', 'dh_fit_dy', 'h_robust_spread', 'h_rms','h_mean','sigma_h_mean','h_med', ...
    'sigma_dh_fit_dx', 'fpb_error' ,'fpb_med_corr','fpb_mean_corr', 'med_r_fit', ...
    'N_initial', 'N_noise', 'n_fit_photons','fpb_N_corr', 'sigma_photon_est', ...
    'TX_med_corr','TX_mean_corr', 'h_LI', 'h_LI_sigma', 'z0','m0', 'x_RGT','y_RGT', ...
    'lat_ctr','lon_ctr', ...
    'fpb_med_corr_sigma','fpb_mean_corr_sigma', ...
    'first_seg_pulse','N_seg_pulses','time','SNR', 'SNR_significance','exit_iteration', ...
    'surface_detection_failed','surface_height_invalid','cross_track_slope_invalid' ,'fpb_correction_failed', ...
    'lat_mean', 'lon_mean', 'delta_x_ATC_mean'};

for kf=1:length(D3_fields)
    D3_empty.(D3_fields{kf})=NaN;
end
int_fields={'N_initial','N_noise','n_fit_photons','N_seg_pulses'};
for k=1:length(int_fields)
    D3_empty.(int_fields{k})=0;
end

% NEW 10/2016: initialize ATL06 status and signal_selection_flags
for field={'signal_selection_source', 'signal_selection_status_confident','signal_selection_status_all', 'signal_selection_status_backup' ...
        'surface_detection_failed','surface_height_invalid','cross_track_slope_invalid' ,'fpb_correction_failed'  }
    D3_empty.(field{1})=uint8(0);
end

D3=repmat(D3_empty, [length(seg_list),2]);

% new 4/2016: setup the output version of the dh_hist.
dh_hist.seg_0=find(mod(seg_list, 10)==0)';
dh_hist.seg_1=min(dh_hist.seg_0+9, length(seg_list));
dh_hist.dh=dh_hist_full.dh;
dh_hist.count=uint16(zeros(numel(dh_hist.dh), numel(dh_hist.seg_0)));
[dh_hist.geo_bin_list, dh_hist.segment_id_list]=deal(NaN(10, numel(dh_hist.seg_0)));
[dh_hist.N_pulses, dh_hist.x_RGT_mean, dh_hist.x_PS_mean, dh_hist.y_PS_mean, dh_hist.lat_mean, dh_hist.lon_mean, dh_hist.bckgrd_per_bin]=deal(NaN(1, numel(dh_hist.seg_0)));
dh_hist=repmat(dh_hist, 1, 2);

tic
for k0=1:length(seg_list)
    % initial processing: find the intial vertical bin centers
    ybar=NaN(1,2);
    clear D2sub
    for kB=1:2
        D3(k0, kB).seg_count=seg_list(k0);
        [D3(k0, kB).x_RGT, x0]=deal(full(dist_for_segment{kB}(seg_list(k0))));
        if D3(k0, kB).x_RGT==0; D3(k0, kB).x_RGT=NaN; continue; end
        AT_els=D2(kB).segment_id==seg_list(k0)-1 | D2(kB).segment_id==seg_list(k0);
        % N.B.  This should be changed to make the N_seg_pulses equal to
        % velocity_sc scaled to ground speed times 
        %D3(k0, kB).N_seg_pulses=57;
         
        
        if isfield(params(kB),'skip') &&  params(kB).skip
            AT_els=[];
        end
        
        if SAVELOG; LOG(k0, kB).AT_els=AT_els; end
        
        if ~any(AT_els)
            D3(k0, kB).signal_selection_source=3;
            D3(k0, kB).signal_selection_status_confident=3;
            D3(k0, kB).signal_selection_status_all=3;
            D3(k0, kB).signal_selection_status_backup=3;
            D3(k0, kB).N_seg_pulses=NaN;
            continue;
        end
        D2sub(kB)=index_struct(D2(kB), AT_els);
        %D3(k0, kB).N_seg_pulses=max(56, diff(range(D2sub(kB).pulse_num))+1);
        %!!! changed this 4/24/19 to match ASAS
        D3(k0, kB).N_seg_pulses=diff(range(D2sub(kB).pulse_num));
        % try choosing the ground strategy
        [initial_fit_els, D3(k0, kB).signal_selection_source, D3(k0, kB).signal_selection_status_confident,  D3(k0, kB).signal_selection_status_all]=...
            choose_ground_strategy(D2sub(kB));
 
        if D3(k0, kB).signal_selection_source < 1  % enough confident PE found to determine slope and height of the bin
            initial_Hwin_min=3;
        else %
            initial_Hwin_min=10;
        end
        
        if D3(k0, kB).signal_selection_source <=1  % if we are using ATL03 flagged PE
            [initial_fit_els, D3(k0, kB).w_surface_window_initial, D3(k0, kB).h_initial, ~]=...
                initial_at_fit(D2sub(kB).x_RGT, D2sub(kB).h_ph, initial_fit_els, x0, median(D2sub(kB).BGR), initial_Hwin_min, D3(k0, kB).N_seg_pulses, params(kB));
        end
        
        % define h_range_initial and (for cases 2 and 3) w_surface_window_initial
        h_range_initial=[NaN NaN];
        if D3(k0, kB).signal_selection_source ==2 &&  D3(k0, kB).signal_selection_status_backup ~=3  % not enough confident or padded PE have been found, fall back to alternate strategies
            [initial_fit_els, D3(k0, kB).signal_selection_source, D3(k0, kB).signal_selection_status_backup]=...
                backup_signal_finding_strategy(D2sub(kB), D2(kB), seg_list(k0), 10);
            h_range_initial(kB)=diff(range(D2sub(kB).h_ph));
            if any(initial_fit_els)
                D3(k0, kB).w_surface_window_initial=diff(range(D2sub(kB).h_ph(initial_fit_els)));  
                D3(k0, kB).h_initial=mean(D2sub(kB).h_ph(initial_fit_els));
            else
                D3(k0, kB).w_surface_window_initial=NaN;
                D3(k0, kB).h_initial=NaN;
            end
        else
            h_range_initial(kB)=D3(k0, kB).w_surface_window_initial;
        end
        if SAVELOG; LOG(k0, kB).initial_fit_els=initial_fit_els; end   
        
        if D3(k0, kB).signal_selection_source<3 % go ahead if we have initial PE
            LS_fit_options=struct( 'Nsigma', 3, 'Hwin_min', 3, 'restrict_fit_to_initial_els', false,'N_it', 25);
            if SAVELOG
                LS_fit_options.SAVELOG=true;
            end
            
            if sum(initial_fit_els) < 10; continue; end
            
            % restrict the fitting to the initial els
            if D3(k0, kB).signal_selection_source<=1
                LS_fit_options.restrict_fit_to_initial_els=true;
            end
            
            if SAVELOG
                [D3(k0, kB), r, els, LOG(k0, kB).LS_fit]=ATLAS_LS_fit(D2sub(kB), D3(k0, kB).x_RGT, initial_fit_els, D3(k0,kB).w_surface_window_initial, params(kB), D3(k0, kB), LS_fit_options);
            else
                [D3(k0, kB), r, els]=ATLAS_LS_fit(D2sub(kB), D3(k0, kB).x_RGT, initial_fit_els, D3(k0,kB).w_surface_window_initial, params(kB), D3(k0, kB), LS_fit_options);
            end
            
            D3(k0, kB).delta_x_ATC_mean=mean(D2sub(kB).x_RGT(els));
            
            if SAVELOG
                LOG(k0, kB).els_after_LS_fit=els;
                LOG(k0, kB).r_after_LS_fit=r;
            end
            
            % calculate the SNR significance
            if exist('SNR_F_table','var')
                D3(k0, kB).SNR_significance=...
                    interpn(SNR_F_table.SNR, SNR_F_table.W_surface_window_initial, SNR_F_table.BGR, SNR_F_table.P_NoiseOnly, ...
                max(min(SNR_F_table.SNR), min(max(SNR_F_table.SNR), D3(k0, kB).SNR)), ...
                max(3, min(max(SNR_F_table.W_surface_window_initial), h_range_initial(kB))), ...
                max(min(SNR_F_table.BGR), min(max(SNR_F_table.BGR) , D3(k0, kB).BGR)));
            end
            
             %NEW: bail out if the SNR_significance is too high
             if D3(k0, kB).SNR_significance >0.05
                 continue
             end
            
            D3(k0, kB).med_r_fit=median(r);
            D3(k0, kB).n_fit_photons=length(r);
            
            % make the time histogram:
            clear temp;
            temp.time=-2/CONSTANTS.c*r;
            %t_WF=((min(temp.time)-2.5*CONSTANTS.dt_hist):CONSTANTS.dt_hist:(max(temp.time)+params(kB).t_dead))';
            % CHANGED 4/2019 to match ASAS
            t_WF=(min(temp.time):CONSTANTS.dt_hist:(max(temp.time)+params(kB).t_dead))';
            N_WF=my_histc(temp.time, t_WF);
            
            % quit if we've somehow lost photons
            if sum(N_WF) < 10 || ~any(isfinite(N_WF)); continue; end
            
            if ~isfield(params,'N_channels'); [params(:).N_channels]=deal(params(:).N_det); end
            % Changed 3/24/2016
            %[D3(k0, kB).fpb_med_corr, D3(k0, kB).fpb_mean_corr, D3(k0, kB).fpb_N_corr, N_WF_corr, D3(k0, kB).fpb_med_corr_sigma, D3(k0, kB).fpb_mean_corr_sigma, minGain, fpb_gain]=...
            %    fpb_corr(t_WF, N_WF, params(kB).N_channels, D3(k0, kB).N_seg_pulses, params(kB).t_dead,  0.01, CONSTANTS.dt_hist);
            [D3(k0, kB).fpb_med_corr, D3(k0, kB).fpb_mean_corr, D3(k0, kB).fpb_N_corr, N_WF_corr, D3(k0, kB).fpb_med_corr_sigma, D3(k0, kB).fpb_mean_corr_sigma, minGain, fpb_gain]=...
                fpb_corr_nonparalyzable(t_WF, N_WF, params(kB).N_channels, D3(k0, kB).N_seg_pulses, params(kB).t_dead, CONSTANTS.dt_hist);
            
            % check if gain correction is valid
            if minGain < 1/(2*params(kB).N_channels)
                D3(k0, kB).fpb_correction_failed=true;
            end
            if D3(k0, kB).signal_selection_source>=2
                D3(k0, kB).surface_detection_failed=1;
            end
            
            D3(k0, kB).N_noise=median(D2sub(kB).BGR)*D3(k0, kB).w_surface_window_final/(CONSTANTS.c/2)*D3(k0, kB).N_seg_pulses;

            if isfield(params(kB),'WF')
                [D3(k0, kB).TX_med_corr, D3(k0, kB).TX_mean_corr, syn_WF]=correct_for_TX_shape(params(kB).WF.t, params(kB).WF.p, D3(k0, kB).w_surface_window_final/(1.5e8),...
                    D3(k0, kB).h_robust_spread/1.5e8, D3(k0, kB).SNR);
            else
                D3(k0, kB).TX_med_corr=0;
            end
            
            if SAVELOG
               LOG(k0, kB).t_WF=t_WF;
               LOG(k0, kB).N_WF=N_WF;
               LOG(k0, kB).N_WF_corr=N_WF_corr;
               LOG(k0, kB).gain=fpb_gain;
               if exist('syn_WF','var');
                   LOG(k0, kB).syn_WF_T=syn_WF.t;
                   LOG(k0, kB).syn_WF_P=syn_WF.P;
                   LOG(k0, kB).syn_WF_mask=syn_WF.mask;
               end
            end
            
            % Calculate the final h_LI
            D3(k0, kB).h_LI=D3(k0, kB).h_mean + D3(k0, kB).fpb_med_corr + D3(k0, kB).TX_med_corr;
  
            if false
                % extra code to estimate the true fpb error
                GG=[ones(size(D2sub_all(kB).h_ph)), real(D2sub_all(kB).x_RGT)-D3(k0, kB).x_RGT];
                els_all=abs(D2sub_all(kB).h_ph-GG*[D3(k0, kB).h_mean; D3(k0, kB).dh_fit_dx]) < D3(k0, kB).w_surface_window_final/2;
                m_all=GG(els_all,:)\D2sub_all(kB).h_ph(els_all);
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
            % pull out the beam number
            temp=unique(D2sub(kB).beam(isfinite(D2sub(kB).beam)));
            D3(k0, kB).beam=temp(1);
            if isfield(D2sub(kB),'track')
                temp=unique(D2sub(kB).track(isfinite(D2sub(kB).track)));
                D3(k0, kB).track=temp(1);
            end
            D3(k0, kB).first_seg_pulse=min(D2sub(kB).pulse_num);
            if isfield(D2sub,'time')
                D3(k0, kB).time=median(D2sub(kB).time);
            else
                D3(k0, kB).time=median(D2sub(kB).delta_time);
            end
            % NEW 10/2015:  Use regression to calculate reference center parameters
            % regress WRT along-track distance to get parameters for
            % segment center: 
            %'x_PS_ctr', 'y_PS_ctr', 'lat_ctr','lon_ctr'
            Ginv_AT=[ones(size(D2sub(kB).x_RGT(els))) D2sub(kB).x_RGT(els)-D3(k0, kB).x_RGT]\eye(sum(els));
            Ginv_AT=Ginv_AT(1,:);
            if isfield(D2sub,'lat')
                D3(k0, kB).lat_ctr=Ginv_AT*D2sub(kB).lat(els);
                D3(k0, kB).lon_ctr=Ginv_AT*D2sub(kB).lon(els);
                %DEBUGGING 9/18
                D3(k0, kB).lat_mean=mean(D2sub(kB).lat(els));
                D3(k0, kB).lon_mean=mean(D2sub(kB).lon(els));
            else
                D3(k0, kB).lat_ctr=Ginv_AT*D2sub(kB).lat_ph(els);
                D3(k0, kB).lon_ctr=Ginv_AT*D2sub(kB).lon_ph(els);
            end
                                 
            if isnan(D3(k0, kB).h_LI)
                D3(k0, kB).surface_height_invalid=true;
            end
            if isfield(D2sub,'y_RGT')
                ybar(kB)=median(D2sub(kB).y_RGT);
            else
                ybar(kB)=median(D2sub(kB).dist_ph_across);
            end
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
        %D3(k0, kB).h_expected_rms=sqrt((dh_fit_dy.^2+D3(k0, kB).dh_fit_dx.^2)*params(kB).sigma_x.^2 + (params(kB).sigma_pulse*1.5e8).^2);
        % 10/2018: CHANGED TO:
        D3(k0, kB).h_expected_rms=sqrt((D3(k0, kB).dh_fit_dx*params(kB).sigma_x).^2 + (params(kB).sigma_pulse*1.5e8).^2);
        D3(k0, kB).sigma_photon_est=sqrt((D3(k0, kB).N_noise*(D3(k0, kB).w_surface_window_final*.287).^2+N_sig*D3(k0, kB).h_expected_rms.^2)/D3(k0, kB).n_fit_photons);
        % CHANGED THIS FROM h_robust_spread to h_rms
        %sigma_per_photon=max(D3(k0, kB).sigma_photon_est, D3(k0, kB).h_robust_spread); 
        sigma_per_photon=max(D3(k0, kB).sigma_photon_est, D3(k0, kB).h_rms);
        D3(k0, kB).sigma_h_mean=D3(k0, kB).sigma_h_mean*sigma_per_photon;
        D3(k0, kB).sigma_dh_fit_dx=D3(k0, kB).sigma_dh_fit_dx*sigma_per_photon;
        D3(k0, kB).h_LI_sigma=max(D3(k0, kB).sigma_h_mean, D3(k0, kB).fpb_med_corr_sigma);
    end
    
    for kB=1:2
        if isfield(params,'RGT')
            D3(k0, kB).RGT=params(kB).RGT;
            D3(k0, kB).PT=params(kB).PT;
            D3(k0, kB).cycle=params(kB).cycle;
        end
    end
    % new, 4/2016: calculate residual histogram
    for kB=1:2
        if ~isfinite(D3(k0, kB).h_mean); continue; end
        xx=D2sub(kB).x_RGT-D3(k0, kB).x_RGT;
        zz=D2sub(kB).h_ph;
        zz=zz(abs(xx)<10);  % select elements from the non-overlapping parts of the segment
        xx=xx(abs(xx)<10);
        if isempty(xx); continue; end
        G=[ones(size(xx)), xx];
        r=zz-G*[D3(k0, kB).h_mean; D3(k0, kB).dh_fit_dx];
        dh_hist_full(kB).count(:, k0)=uint8(my_histc(r, dh_hist_full(kB).dh));
    end
    
    if SAVELOG
        for kB=1:2
            LOG(k0, kB).D3=D3(k0, kB);
            if exist('D2sub','var')
                LOG(k0, kB).D2=D2sub(kB);
            end
        end
    end
    
    
%     if mod(k0, 250)==0
%         disp([num2str(k0) ' out of ' num2str(length(seg_list))]);
%     end
end

if isempty(D3)
    D3a=[];
    return
end

clear D3a;
f=fieldnames(D3); for kF=1:length(f); for kB=1:2; D3a.(f{kF})(:,kB)=cat(1, D3(:, kB).(f{kF})); end; end



%  new 4/2016 : collect the D3 to D3a:
clear D3a;
ff=fieldnames(D3);
for kF=1:length(ff)
    for kB=1:2
        D3a.(ff{kF})(:,kB)=cat(1, D3(:, kB).(ff{kF}));
    end
end

D3a.ATL06_quality_summary=calc_ATL06_summary_flag(D3a);


% new 4/2016: stack the dh_histogram by 10.
for kB=1:2
    for ks=1:length(dh_hist(kB).seg_0)
        els=dh_hist(kB).seg_0(ks):dh_hist(kB).seg_1(ks);
        
        %els=els(D3a.surface_height_invalid(els, kB)==0   & D3a.SNR(els, kB) > 0.5);
        els=els(D3a.ATL06_quality_summary(els, kB)==0);
        if ~isempty(els)
            dh_hist(kB).bckgrd_per_bin(ks)=sum(D3a.BGR(els, kB).*D3a.N_seg_pulses(els, kB)/2*CONSTANTS.residual_histogram_spacing/(CONSTANTS.c/2));
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
function [initial_els, w_surface_window, h_window_ctr, AT_slope]=initial_at_fit(x_AT, h, initial_els, x0, BGR, W_min, N_pulses_in_seg, params)
c2=3e8/2;

G=[ ones(size(x_AT)), x_AT-x0];
m=G(initial_els,:)\h(initial_els);
h_window_ctr=m(1); AT_slope=m(2);
r=h(:)-G*m;
Noise_Ph_per_m=N_pulses_in_seg*BGR/c2;
 
H_win=diff(range(r(initial_els)));
sigma_r=robust_peak_width_CDF(r(initial_els), Noise_Ph_per_m*H_win, [0 H_win]-H_win/2);
sigma_expected=sqrt((params.sigma_x*m(2)).^2+ (params.sigma_pulse*c2).^2);
w_surface_window=max(W_min, 6*max(sigma_r, sigma_expected));
initial_els= abs(r) < w_surface_window/2;

%---------------------------------------------------------------------------------------
function [D3, r0, els, LOG]=ATLAS_LS_fit(D2, L0, initial_els,  H_win, params, D3, options)

if isfield(options,'SAVELOG')
    SAVELOG=true;
else
    SAVELOG=false;
end
if ~exist('options','var')
    options=struct('Nsigma', 3, 'Hwin_min', 3);
end
if ~isfield(options,'N_it')
    options.N_it=25;
end

if options.restrict_fit_to_initial_els
    selectable_els=initial_els;
else
    selectable_els=true(size(D2.h_ph));
end


c2=3e8/2;

if SAVELOG
    LOG.iterations=struct('els', [], 'sigma_r', [], 'sigma_expected', [], 'W_win', [], 'h_ctr', [], 'dhdx', []);
end

% fit an along-track polynomial
%s_ctr=mean(range(real(D2.x_LC)));
ds=real(D2.x_RGT)-L0;
els=initial_els;

G=[ones(size(ds(:))), ds(:)];
m=G(els,:)\D2.h_ph(els);

D3.N_initial=sum(els);
if SAVELOG; LOG.G=G; end
% iterate to reduce residuals
Noise_Ph_per_m=D3.N_seg_pulses*median(D2.BGR)/c2;
if ~exist('H_win','var')|| isempty(H_win)
    H_win=diff(range(D2.h_ph(els)));
end

filter_hist=false(options.N_it, numel(D2.h_ph));

for k=1:options.N_it
    m_last=m;
    m=G(els,:)\D2.h_ph(els);
    
    r_all=D2.h_ph-G*m;
    r0=r_all(els);
    sigma_r=min(5,robust_peak_width_CDF(r0, Noise_Ph_per_m*H_win, [0 H_win]-H_win/2));
    sigma_expected=sqrt((c2*params.sigma_pulse).^2+params.sigma_x.^2*(m(2).^2));
    
    els_last=els;
    filter_hist(k,:)=els;
    SNR_last=sum(els)/(H_win*Noise_Ph_per_m);
    H_win_last=H_win;
    
    H_win=max([2*[sigma_expected, sigma_r]*options.Nsigma, 0.75*H_win_last, options.Hwin_min]);
    els=selectable_els & abs(r_all ) < H_win/2;
    if SAVELOG
        LOG.iterations(k).els=els;
        LOG.iterations(k).sigma_r=sigma_r;
        LOG.iterations(k).sigma_expected=sigma_expected;
        LOG.iterations(k).W_win=H_win;
        LOG.iterations(k).h_ctr=m(1);
        LOG.iterations(k).dhdx=m(2);
    end
     
    if sum(els) < 10 || diff(range(D2.x_RGT(els))) < 20   
        m=m_last;
        H_win=H_win_last;
        els=els_last;
        break
    end
    
    if sum(abs(els_last-els))==0  % no change from last iteration, or we've converged.
        break
    end
 
end


D3.h_mean=m(1);
D3.dh_fit_dx=m(2);
D3.BGR=median(D2.BGR);
r0=r_all(els);
D3.w_surface_window_final=H_win;

D3.h_robust_spread=iqr(r0)/2;
D3.h_rms=std(r0);
D3.h_med=m(1)+median(r0);

% NEW 10/2015: Calculate the noise estimate and the SNR
Noise_est=H_win*Noise_Ph_per_m;
D3.SNR=(sum(els)-Noise_est)./Noise_est;

% NEW 10/2015: Report the exit iteration
D3.exit_iteration=k;
if options.N_it==1
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
els=D2.signal_conf_ph > 1 &  isfinite(D2.h_ph) ;

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
els=D2.signal_conf_ph >= 1 &  isfinite(D2.h_ph);
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
function [selected_PE, signal_selection_source, signal_selection_status]=backup_signal_finding_strategy(D2_local, D2_all, segment_id, W)

found=D2_local.signal_conf_ph > 1 &   isfinite(D2_local.h_ph);
if any(found)
    % if we have any selected PE, center the window on them, see if this
    % gives us a valid window (of any quality
    Hmean=mean(D2_local.h_ph(found));
    els=abs(D2_local.h_ph-Hmean)<5;
    if sum(els) > 10 && diff(range(D2_local.x_RGT(els)))>20
        selected_PE=els;
        signal_selection_source=2;
        signal_selection_status=0;
        return
    end
end

% nothing found, or the PE around the found PE were not usable.
these=D2_all.segment_id >=  segment_id-2 & D2_all.segment_id <= segment_id+1;
if sum(these)<10
    selected_PE=[];
    signal_selection_source=3;
    signal_selection_status=4;
    return
end

% use the histogram strategy
h=D2_all.h_ph(these);
bins=(floor(min(h))+0.25):0.5:ceil(max(h));

C1=zeros(size(bins));
for kB=1:length(C1); 
    C1(kB)=sum(abs(h-bins(kB))<W/2);
end

% procede if more than 20 PE are in the best bin
if max(C1) > 20   
    % select the bins that are not significantly differnt from the peak
    z0r=range(bins(C1>max(C1)-sqrt(max(C1))));
    selected_PE=D2_local.h_ph >= z0r(1)-W/2 & D2_local.h_ph <= z0r(2)+W/2;
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

%---------------------------------------------------------------
function flag=calc_ATL06_summary_flag(D3)

flag=~(D3.h_robust_spread<1 &...
    D3.signal_selection_source <=1 & ...
    D3.h_LI_sigma < 1 & ...
    D3.SNR_significance < 0.02);



