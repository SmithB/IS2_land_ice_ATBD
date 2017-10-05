function varargout=ATL11_proc_ATBD_Nseg(varargin)

if nargout>0
    [varargout{1:nargout}]=feval( varargin{:});
else
    feval(varargin{:});
end
%------------------------------------------
function [D3b, params]=calc_fit(D3a, N_reps,  params)

if isfield(params,'DOPLOT')
    if ~isstruct(params.DOPLOT)
        params.DOPLOT={'xtrack_surf','input_data_edit', 'h_of_t'};
    end
    clf
    params.axes=cheek_by_jowl(2, length(params.DOPLOT), [0 0 1 1]*.6+0.2*[1 1 0 0]);
end

if ~isfield(params,'N_segs')
    params.N_segs=5;
end

params.AT_search=20*params.N_segs/2+2;

% set defaults
params.beam_spacing=90;
params.y_scale=100;
params.ref_seg_Lsearch=125;

% work out how the polynomial components will be reported:
max_deg=struct('x', params.max_degree_x, 'y', params.max_degree_y);
[~, deg_ind]=sort_poly2_components([], [],  max_deg);
params.poly_exp_x=deg_ind(:, 1); 
params.poly_exp_y=deg_ind(:, 2);

%corrected_h subgroup
fields.corrected_h={'ref_pt_lat',  'ref_pt_lon', 'mean_pass_time', ...
    'pass_h_shapecorr', 'pass_h_shapecorr_sigma','pass_h_shapecorr_sigma_systematic', ...
    'pass_included_in_fit','pass_quality_summary'};
% reference_points group
fields.reference_point={'PT', 'RGT', 'L_search', 'rgt_azimuth', 'x_ATC', 'y_ATC'};
% reference surface group
fields.ref_surf={  'slope_change_rate','slope_change_rate_sigma','poly_ref_surf', 'poly_ref_surf_sigma', ...
    'surf_fit_misfit_RMS','surf_fit_misfit_chi2', ...
    'N_cycle_avail', 'N_cycle_fit', ...
    'fit_X_slope', 'fit_Y_slope','fit_curvature', 'ATL11_status', 'P_value'};
% pass quality stats
fields.pass_quality={'ATL06_summary_zero_count','min_SNR_significance','mean_uncorr_reflectance','min_signal_selection_source'};
% pass stats
fields.pass_stats={'h_robust_spread_mean', 'h_rms_mean', 'reflct_uncorr_mean', 'tide_ocean_mean','mean_pass_lat',  'mean_pass_lon', ...
    'h_robust_spread_mean', 'cloud_flg_best', 'bsl_h_mean', 'bslh_conf_best', 'laser_beam', 'y_RGT', 'ref_pt_number', ...
    'sigma_g_h', 'sigma_g_x', 'sigma_g_y', 'snr_mean'};
% fields that aren't on ATL11 but are useful here 
fields.non_product={'z0' 'delta_z_signal', 'delta_z_error','delta_z_firn','x_PS_ctr', 'y_PS_ctr'};
 
f0=fieldnames(fields);
for kF0=1:length(f0)
    clear temp;
    for kF1=1:length(fields.(f0{kF0}))
        temp.(fields.(f0{kF0}){kF1})=[];
    end
    D3b.(f0{kF0})=temp;
end

ATL11_status=uint8(0);

valid=input_data_edit(D3a, params);
if ~any(valid.combined); return;  end

fit_rep_x0_list=unique(valid.Rep_x0(valid.combined,:),'rows');
if length(unique(fit_rep_x0_list(:,1))) < 2
    return
end
Rep_x0=[D3a.rep(:), D3a.x_RGT(:)];
fit_segs=ismember(Rep_x0, fit_rep_x0_list, 'rows');
% 
% % copy the segment stats
% for field=fields.segment_stats
%     if isfield(D3a, field{1});
%         D3b.segment_stats.(field{1})=D3a.(field{1})(:);
%     end
% end
% D3b.segment_stats.seg_count=zeros(size(D3a.h_LI))+params.ref_pt_num;
% 
% % copy the segment geolocation
% for field=fields.segment_geolocation
%     if isfield(D3a, field{1})
%         D3b.segment_geolocation.(field{1})=D3a.(field{1})(:);
%     end
% end
for field=fields.ref_surf
    if isfield(D3a, field{1})
        D3b.ref_surf.(field{1})=D3a.(field{1})(:);
    end
end

% pick the degree:
DOF_y=length(unique(fit_rep_x0_list(:,1)))+1;
poly_degree.y=min(params.max_degree_y, max(1, DOF_y-2));
DOF_x=length(unique(fit_rep_x0_list(:,2)));
poly_degree.x=min(params.max_degree_x, max(1, DOF_x-2));

% pick the y center
[fit_ctr.y, valid]=choose_y0(D3a, valid, params);

% pick the x center.  If the data span the segment center, use it
if min(D3a.x_RGT(fit_segs))< params.x_RGT_ctr && max(D3a.x_RGT(fit_segs)) > params.x_RGT_ctr
    fit_ctr.x=params.x_RGT_ctr;
else
    % otherwise, use the mean of the segment locations
    fit_ctr.x=mean(range(D3a.x_RGT(fit_segs)));
end

% if there is a time span longer than 1.5 years, fit for the delta_slope
% for all components for which we're recovering the slope.
if diff(range(D3a.time(fit_segs))) >= 1.5
    if poly_degree.x >=1
        poly_degree.dsxdt=1;
    else
        poly_degree.dsxdt=0;
    end
    if poly_degree.y >=1
        poly_degree.dsydt=1;
    else
        poly_degree.dsydt=0;
    end
else
    poly_degree.dsxdt=0;
    poly_degree.dsydt=0;
end

if sum(valid.combined)==0; ATL11_status=bitor(ATL11_status, 1); end

% fit the data 
[CH, RS,  poly_vals, valid]=est_sigma_and_hrep(D3a, 1:N_reps, valid, fit_ctr, poly_degree, ATL11_status, params);

% assign corrected_h subgroup
D3b.reference_point.x_ATC=fit_ctr.x;
D3b.reference_point.y_ATC=fit_ctr.y;
D3b.reference_point.rgt_azimuth=mean(D3a.azimuth(isfinite(D3a.azimuth)));   % should include a 360-degree wrap here

fit_rep_list=find(CH.pass_included_in_fit); %find(~bitand(uint8(CH.cycle_selection_field), 16));
% put together the mean pass time, etc
[CH.mean_pass_time]=reduce_by(@mean, D3a,  {'time'},'rep', 1:N_reps);
ref_pt_subset=index_struct(D3a, ismember(D3a.rep, fit_rep_list) & isfinite(D3a.x_RGT) & isfinite(D3a.y_RGT),{'lat_ctr','lon_ctr','x_RGT', 'y_RGT'});
[CH.ref_pt_lat, CH.ref_pt_lon]=regress_to(ref_pt_subset, {'lat_ctr', 'lon_ctr'}, {'x_RGT','y_RGT'}, [fit_ctr.x, fit_ctr.y]);
D3b.corrected_h=CH;

% calculate weighted mean, etc fields in pass_stats
weights=D3a.h_LI_sigma.^(-2); 
weights(~valid.selected_segs)=NaN;
D3b.pass_stats=assign_pass_stats( D3a, weights,  1:N_reps, fields.pass_stats);

% assign the pass quality stats
[D3b.pass_quality, pass_quality_summary]=assign_pass_quality_stats(N_reps, D3a, valid, fields.pass_quality );
D3b.corrected_h.pass_quality_summary=pass_quality_summary;

% assign the non-product fields
[D3a.x_PS_ctr, D3a.y_PS_ctr]=ll2ps(D3a.lat_ctr, D3a.lon_ctr);
for k=1:length(fields.non_product)
    D3b.non_product.(fields.non_product{k})= reduce_by(@mean, D3a,  fields.non_product(k) ,'rep', 1:N_reps);
end

% assign reference point group outside of this function
ff=fieldnames(RS);
for kf=1:length(ff)
    D3b.ref_surf.(ff{kf})=RS.(ff{kf});
end
% plot the plots that need the plotting
if isfield(params,'DOPLOT') && ismember('xtrack_surf', params.DOPLOT );
    xtrack_plot(D3a, hc_rep, y0, mean(hc_rep(ismember(1:N_reps, fit_rep_list))), m_poly, sigma_m_poly, fit_rep_list, params)
    if size(params.axes,2)>1;
        set(params.axes(:,end),'yaxislocation','right');
    end
    reps=1:N_reps;
    fprintf(1, 'rep | mx | my | y | A\n');
    fprintf(1, '_______________________\n');
    fprintf('%3d |  %1d   %1d    %1d   %1d\n', [reps; ismember(reps, good_x_reps); ismember(reps, good_y_reps); ismember(reps, good_dist_reps); ismember(reps, fit_rep_list)]);
end

if isfield(params,'DOPLOT') && ismember('h_of_t', params.DOPLOT);
    axes(params.axes(1,end));
    P1=get(params.axes(1, end),'position');
    P2=get(params.axes(2,end), 'position');
    delete(params.axes(1, end)); 
    set(params.axes(2,end),'position', [P1(1)+.1 P2(2), P1(3), P1(2)+P1(4)-P2(2)]);
    axes(params.axes(2,end));
    plot(D3b.corrected_h.mean_pass_time/365, D3b.corrected_h.pass_h_shapecorr,'r.');
    errorbar(D3b.corrected_h.mean_pass_time/365, D3b.corrected_h.pass_h_shapecorr, D3b.corrected_h.pass_h_shapecorr_sigma,'ko');
    set(gca,'yaxislocation','right');
    ylabel('Corrected h, m');
    set(params.axes(:,2),'yaxislocation','right');
    set(params.axes(1, 1:2),'xtick',[]);
    xlabel(params.axes(2,1),'cross-track y, m');
    xlabel(params.axes(2,2),'cross-track y, m');
end
    


%-------------------------------------------------
function [CH, RS, SV, valid, params_out]=est_sigma_and_hrep(D3a, reps, valid,  fit_ctr, poly_degree, ATL11_status, params)
 
% calculate the polynomial correction for a set of segments.
% outputs: 
% CH: Corrected heights for reps
% RS: Reference surface
% SV: segment values
% valid: set of validity flags

deg_ind=[params.poly_exp_x(:), params.poly_exp_y(:)];
%  [hc_rep, sigma_hc_rep, sigma_hc_rep_systematic, m_surf, sigma_m_surf, h_poly_seg, sigma_h_poly_seg, r_fit, chi2_r]
fields.perPass={'mean_pass_time', 'pass_h_shapecorr', 'pass_h_shapecorr_sigma','pass_h_shapecorr_sigma_systematic'};
valid.selected_segs=false(size(D3a.h_LI));

fields.SV={'h_poly_seg','sigma_h_poly_seg','sigma_hc_seg_systematic'};
 
for kf=1:length(fields.perPass)
    CH.(fields.perPass{kf})=NaN(1, length(reps));
end
for kf=1:length(fields.SV)
    SV.(fields.SV{kf})=NaN(size(D3a.time)); 
end
SV.accepted_by_iterative_fit=false(size(D3a.time));

% assign the output values for the poly-fit matrix
[~, poly_coeff_degree]=poly2_fit_mat([], [], struct('x', params.max_degree_x, 'y', params.max_degree_y));
[RS.slope_change_rate, RS.slope_change_rate_sigma]=deal(NaN(1,2));
[RS.poly_ref_surf, RS.poly_ref_surf_sigma]=deal(NaN(1, length(poly_coeff_degree.x)));
[RS.surf_fit_misfit_RMS, RS.surf_fit_misfit_chi2, RS.N_cycle_avail, RS.N_cycle_fit, RS.fit_X_slope, RS.fit_Y_slope, RS.fit_curvature]=deal(NaN);
 
RS.ATL11_status=ATL11_status;

t0=params.time_zero;
fit_Rep_x0=valid.Rep_x0(valid.combined,:);
fit_reps=unique(fit_Rep_x0(:,1));
N_fit_reps=length(unique(fit_Rep_x0(:,1)));

Rep_x0=[D3a.rep(:), D3a.x_RGT(:)];
ref_fit_seg=ismember(Rep_x0, fit_Rep_x0,'rows');

G_full_z0=zeros(size(D3a.h_LI(:),1),N_fit_reps);
for kr=1:length(fit_reps)
    G_full_z0(ismember(Rep_x0, fit_Rep_x0,'rows') & D3a.rep==fit_reps(kr), kr)=1;
end

[G_full_poly, poly_coeff_degree]=poly2_fit_mat((D3a.x_RGT(:)-fit_ctr.x)/params.y_scale, (D3a.y_RGT(:)-fit_ctr.y)/params.y_scale, poly_degree);

% fit for change in slope iff we are fitting for it.
G_full_dsdt=[];
if poly_degree.x > 0
    G_full_dsdt=[G_full_dsdt (D3a.x_RGT(:)-fit_ctr.x)/params.y_scale.*(D3a.time(:)-t0)/params.t_scale];
end
if poly_degree.y > 0
    G_full_dsdt=[G_full_dsdt (D3a.y_RGT(:)-fit_ctr.y)/params.y_scale.*(D3a.time(:)-t0)/params.t_scale];
end

% build the table of contents for G
TOC.cols.dsdt=1:size(G_full_dsdt,2); 
TOC.cols.poly=TOC.cols.dsdt(end)+(1:size(G_full_poly,2));
TOC.cols.poly_degree=poly_coeff_degree;
TOC.cols.z0=TOC.cols.poly(end)+(1:size(G_full_z0,2));

sigma_seg=sqrt(D3a.h_LI_sigma.^2+params.sigma_geo.^2*(D3a.dh_fit_dx.^2+D3a.dh_fit_dy.^2));

G=[G_full_dsdt, G_full_poly, G_full_z0];
this_fit_seg=ref_fit_seg(:);
FailedFit=false;
N_iterations_max=5;
for k=1:N_iterations_max
    % iterate to eliminate outliers
    %subset G
    sz=sum(this_fit_seg);       
    G1=G(this_fit_seg,:);
    % only solve for columns of G that have a non-zero range
    cols=max(G1)-min(G1)>0;
    if sum(cols)==0; RS.ATL11_status=bitor(RS.ATL11_status, 4); FailedFit=true; break;  end
    G1=G1(:,cols);
    % quit if the solution is underdetermined
    if size(G1,2)>size(G1,1); RS.ATL11_status=bitor(RS.ATL11_status, 4); FailedFit=true; break;  end
    Cinv_sub=spdiags(sigma_seg(this_fit_seg).^-2, 0, sz, sz);
    
    %calculate Ginv and the fit
    if condest(G1'*Cinv_sub*G1) > 1e13
        FailedFit=true;
        break
    end
    G1inv=(G1'*Cinv_sub*G1)\G1';
    m1=zeros([size(G,2),1]);
    m1(cols)=G1inv*Cinv_sub*D3a.h_LI(this_fit_seg);
    
    % calculate the misfit statistics:
    r_fit=D3a.h_LI(this_fit_seg)-G1*m1(cols);
    this_sigma=iqr(r_fit./sigma_seg(this_fit_seg))/2;
    RS.surf_fit_misfit_chi2=(r_fit'*Cinv_sub*r_fit);
    % P value--- ~1 means a good misfit
    RS.P_value=1-chi2cdf(RS.surf_fit_misfit_chi2, size(G1,1)-size(G1,2));
    
    if  RS.P_value < 0.025 && k < N_iterations_max
        % iterate if we're too far from a reasonable chi2 value
        this_r=D3a.h_LI-G*m1;
        next_fit_seg = ref_fit_seg & abs(this_r./sigma_seg) < 3*this_sigma;
        if all(next_fit_seg==this_fit_seg); this_fit_seg=next_fit_seg; continue; end
        this_fit_seg=next_fit_seg;
    else
        continue
    end
end

if FailedFit
    r_fit=NaN;
    [sigma_m, m1]=deal(NaN(size(G(1,:)))');
    this_fit_seg=false(size(this_fit_seg));
    Cm=NaN(size(G,2));
end

RS.surf_fit_misfit_RMS=sqrt(mean(r_fit.^2));
% report the seg selection
valid.iterative_edit_seg=this_fit_seg;
valid.iterative_edit_pair=ismember(valid.Rep_x0,[D3a.rep(this_fit_seg), D3a.x_RGT(this_fit_seg)],'rows');
valid.iterative_edit_rep=ismember(valid.Rep_x0(:,1), D3a.rep(this_fit_seg));

% calculate the surface curvature and slope
XR=round_to(range(D3a.x_RGT(:)-fit_ctr.x), 10)+[-10 10];
YR=round_to(range(D3a.y_RGT(:)-fit_ctr.y), 10)+[-10 10];
[xg,yg]=meshgrid(XR(1):10:XR(2), YR(1):10:YR(2)); 
h_samp=poly2_fit_mat(xg/params.y_scale, yg/params.y_scale, poly_degree)*m1(TOC.cols.poly);
[gx, gy]=gradient(reshape(h_samp, size(xg)), xg(1,:), yg(:,1)); 

ctr_els=abs(xg+1i*yg)<50;
RS.fit_X_slope=mean(gx(ctr_els));
RS.fit_Y_slope=mean(gy(ctr_els));
RS.fit_curvature=sqrt(std(gx(ctr_els)).^2+std(gy(ctr_els)).^2); 
seg_y_slope=interp2(xg(1,:), yg(:,1), gy, D3a.x_RGT(:)-fit_ctr.x, D3a.y_RGT(:)-fit_ctr.y);
seg_x_slope=interp2(xg(1,:), yg(:,1), gx, D3a.x_RGT(:)-fit_ctr.x, D3a.y_RGT(:)-fit_ctr.y);

% assign correlated errors in Cd for each pass, propagate the
% errors to Cm
if ~FailedFit
    [seg_r, seg_c]=ndgrid(D3a.rep(this_fit_seg));
    %Cd_diag=spdiags(max(D3a.h_LI_sigma ,RS.surf_fit_misfit_RMS).^2, 0, sum(this_fit_seg), sum(this_fit_seg));
    Cd_diag=spdiags(max(D3a.h_LI_sigma(this_fit_seg) ,RS.surf_fit_misfit_RMS).^2, 0, sum(this_fit_seg), sum(this_fit_seg));
    Cd_off_diag=params.sigma_geo.^2.*(RS.fit_X_slope.^2+RS.fit_Y_slope.^2 + RS.fit_curvature.^2)*(seg_r==seg_c);
    Cd_tot=Cd_diag+Cd_off_diag;
    Cm=G1inv*Cinv_sub*Cd_tot*Cinv_sub*G1inv';
    sigma_m=NaN(size(G(1,:)))';
    sigma_m(cols)=sqrt(diag(Cm));
end

% report the model and its errors:
surf_cols=[TOC.cols.dsdt, TOC.cols.poly];
G_poly_sub=G(:, surf_cols);

% Rescale the polynomial coefficients, assign to fixed columns, report with
% errors
max_deg=struct('x', params.max_degree_x, 'y', params.max_degree_y);
poly_vals=m1(TOC.cols.poly)';
RS.poly_ref_surf=sort_poly2_components(poly_vals, poly_coeff_degree,  max_deg, deg_ind);
poly_sigma=sigma_m(TOC.cols.poly)';
RS.poly_ref_surf_sigma=sort_poly2_components(poly_sigma, poly_coeff_degree, max_deg, deg_ind);

RS.slope_change_rate=m1(TOC.cols.dsdt)'/params.y_scale/params.t_scale;
RS.slope_change_rate_sigma=sigma_m(TOC.cols.dsdt)'/params.y_scale/params.t_scale;;
RS.N_cycle_avail=length(unique(D3a.rep));
RS.N_cycle_fit=length(unique(valid.Rep_x0(valid.combined,1)));

% clip out the polynomial part of Cm
Cm1=Cm(surf_cols, :);
Cm1=Cm1(:,surf_cols);

% calculate the corrected h for all segs
% first calculate the reference surface for all segs
SV.h_poly_seg=G_poly_sub*m1(surf_cols);
 
% propagate the shape-correction errors back to the segment locations
C_surf=G_poly_sub*Cm1*G_poly_sub';
SV.sigma_h_poly_seg=sqrt(diag(C_surf));

% report the error in the hc params specified by the fit
CH.pass_h_shapecorr_sigma(ismember(reps, fit_reps))=sigma_m(TOC.cols.z0);
CH.pass_h_shapecorr(ismember(reps, fit_reps))=m1(TOC.cols.z0);

%calculate the surface y derivative for error propagation
%yfit_m_components=ismember([poly_coeff_degree.x(:), poly_coeff_degree.y(:)], [ypoly_coef_degree.x(:), ypoly_coef_degree.y(:)],'rows');
%dhdy=G_full_poly_derivative_y*RS.poly_ref_surf(yfit_m_components)';
dhdx=D3a.dh_fit_dx;
%SV.sigma_hc_seg_systematic=params.sigma_geo*sqrt(dhdy(:).^2+dhdx(:).^2);
SV.sigma_hc_seg_systematic=params.sigma_geo*sqrt(seg_y_slope.^2+dhdx(:).^2);
for k=1:length(fit_reps)
    these=(D3a.rep==fit_reps(k) & isfinite(SV.sigma_hc_seg_systematic) & valid.selected_segs);
    CH.pass_h_shapecorr_sigma_systematic(fit_reps(k))=sqrt(mean(SV.sigma_hc_seg_systematic(these).^2));
end

rep_fit_mask=ismember(D3a.rep, fit_reps);
% subset the polynomial matrix for non-fit data:
if any(~rep_fit_mask)    
    % calculate corrected h for all segments, and the errors
    hc_seg=D3a.h_LI-SV.h_poly_seg;
    sigma_hc_seg=sqrt(D3a.h_LI_sigma.^2+SV.sigma_h_poly_seg.^2);
    uRep=reps(~ismember(reps, fit_reps));
    % loop over non-fit reps
    for kR=1:length(uRep)
        % find the output column matching this rep
        this_col=find(reps==uRep(kR));  
        % find the segments in this rep
        these_D3_els=find(D3a.rep==uRep(kR));
        these_errors=sigma_hc_seg(these_D3_els);
        
        %If there are multiple elevations for this rep, pick the one with
        %the best elevation
        best=find(these_errors==min(these_errors), 1, 'first');
        if ~isempty(best)
            valid.selected_segs(these_D3_els(best))=true;
            CH.pass_h_shapecorr(this_col)=hc_seg(these_D3_els(best));
            CH.pass_h_shapecorr_sigma(this_col)=these_errors(best);
            CH.pass_h_shapecorr_sigma_systematic(this_col)=SV.sigma_hc_seg_systematic(these_D3_els(best));
        end
    end
end

SV.accepted_by_iterative_fit=valid.iterative_edit_seg;
valid.selected_segs=valid.selected_segs | valid.iterative_edit_seg; 

for k=1:length(reps)
    these=(D3a.rep==reps(k) & isfinite(SV.sigma_hc_seg_systematic) & valid.selected_segs);
    if any(these)
        CH.pass_h_shapecorr_sigma_systematic(reps(k))=sqrt(mean(SV.sigma_hc_seg_systematic(these).^2));
    end
end

CH.pass_included_in_fit=false(size(reps));
CH.pass_included_in_fit(ismember(reps, valid.Rep_x0(valid.combined & valid.iterative_edit_rep)))=true;


% 
% % calculate the cycle_selection_field bits
% cycle_selection_field=uint8(zeros(size(reps)));
% these=~ismember(reps, D3a.rep); cycle_selection_field(these)=bitor(cycle_selection_field(these), 1);
% these=~ismember(reps, valid.Rep_x0(valid.slope_x)); cycle_selection_field(these)=bitor(cycle_selection_field(these), 2);
% these=~ismember(reps, valid.Rep_x0(valid.slope_y)); cycle_selection_field(these)=bitor(cycle_selection_field(these), 4);
% these=~ismember(reps, valid.Rep_x0(valid.dist)); cycle_selection_field(these)=bitor(cycle_selection_field(these), 8);
% these=~ismember(reps, valid.Rep_x0(valid.combined) ); cycle_selection_field(these)=bitor(cycle_selection_field(these), 16);
% these=~ismember(reps, valid.Rep_x0(valid.iterative_edit_rep) ); cycle_selection_field(these)=bitor(cycle_selection_field(these), 32);
%CH.cycle_selection_field=cycle_selection_field;




%-------------------------------------------------
function [y0, valid]=choose_y0(D3a, valid, params)

% select the y0 value that allows the maximum number of valid repeats 
%to be used in the across-track fit, select all repeats within
%ref_seg_Lsearch of the center.
%
% algorithm: for each repeat, calculate a y range.
% for a range of y0_ref values between min(y_range)+beam_spacing/2 and
% max(y_range-beam_spacing/2, count the number of repeats that fall entirely
% within y0_ref-ref_seg_Lsearch/2 and y0_ref+ref_seg_Lsearch/2, pick the
% maximum
% if there is a tie: break the tie by picking the median y0_ref value
%
% required inputs :
% structure D3a with fields rep, x_RGT, y_RGT
%
% structure valid with fields x0, combined, slope_x
%
% structure params with fields ref_seg_Lsearch, beam_spacing

Rep_x0= [D3a.rep(:), D3a.x_RGT(:)];
fit_rep_x0_list=valid.Rep_x0(valid.combined,:);

% choose y0
fit_seg=ismember(Rep_x0, fit_rep_x0_list,'rows');
y_vals=min(D3a.y_RGT(fit_seg)):2:max(D3a.y_RGT(fit_seg));
rep_count=zeros(size(y_vals));
% 
% uReps=unique(fit_rep_x0_list(:,1));
% for kR=1:length(uReps);
%     YR=range(D3a.y_RGT(D3a.rep==uReps(kR)));
%     these=y_vals >=YR(1) & y_vals <= YR(2);
%     rep_count(these)=rep_count(these)+1;
% end
% y0=median(y_vals(rep_count==max(rep_count)));
for k=1:size(fit_rep_x0_list,1)
    these=D3a.rep==fit_rep_x0_list(k);
    if diff(range(D3a.beam(these)))>0
        this_y=D3a.y_RGT(these);
        YR(k,:)=range(this_y(:));
    end
end
YR=YR(YR(:,1)~=0,:);

for kY=1:length(y_vals)
    rep_count(kY)=sum(sum(abs(YR-y_vals(kY))<=params.ref_seg_Lsearch,2)==2);
end
y0=median(y_vals(rep_count==max(rep_count)));

valid.dist=false(size(valid.slope_x));
for k=1:size(fit_rep_x0_list,1)
    these=D3a.rep==fit_rep_x0_list(k);
    this_y=D3a.y_RGT(these);
    
    these(these)=abs(this_y-y0)<=params.ref_seg_Lsearch;
    if any(these(:)) && any(D3a.beam(these)==1) && any(D3a.beam(these)==2)
        valid.dist(k)=true;
    end
    
%     if min(abs(D3a.y_RGT(D3a.rep==fit_rep_x0_list(k))-y0))<params.ref_seg_Lsearch-params.beam_spacing/2;
%         valid.dist(k)=true;
%     end
end

%-------------------------------------------------
function [valid]=input_data_edit(D3a, params)

if isfield(params,'DOPLOT') && ismember('input_data_edit',params.DOPLOT  )
    DOPLOT=true;
else
    DOPLOT=false; 
end

% data-based edits 
Vmat=[D3a.n_fit_photons >=10  D3a.signal_selection_source <=2  D3a.ATL06_status==0  isfinite(D3a.h_LI)  isfinite(D3a.h_LI_sigma) D3a.SNR_significance < 0.05];
Vmat=[Vmat D3a.h_LI_sigma < 3*max(0.2,median(D3a.h_LI_sigma(all(Vmat,2))))];
%seg_valid.data=D3a.n_fit_photons >=10 & D3a.signal_selection_flags < 32 & D3a.ATL06_status==0 & isfinite(D3a.h_LI) & isfinite(D3a.h_LI_sigma);
%seg_valid.data=seg_valid.data & D3a.h_LI_sigma < 3*median(D3a.h_LI_sigma(seg_valid.data));
seg_valid.data=all(Vmat,2);
y_scale=100;

Rep_x0=[D3a.rep, D3a.x_RGT];
uRep_x0=unique(Rep_x0,'rows');

% loop over ATL06 seg pairs, count valid data segs, select representative x
% and y values
[dhdy, ybar,   sigma_dhdy, xpair, r_slope_y]=deal(NaN(size(uRep_x0,1),1));
[valid.data, valid.slope_x, valid.slope_y]=deal(false(size(uRep_x0,1),1));
for kR=1:size(uRep_x0,1)
    these=find(D3a.rep==uRep_x0(kR,1) & D3a.x_RGT==uRep_x0(kR,2));
    % check if we have a pair
    if length(these)==2
        dhdy(kR)=mean(D3a.dh_fit_dy(these));
        ybar(kR)=mean(D3a.y_RGT(these));
        sigma_dhdy(kR)=median(D3a.sigma_dh_fit_dy(these));
        xpair(kR)=mean(D3a.x_RGT(these));
        if sum(seg_valid.data(these))==2
            valid.data(kR)=true;
        end
    end
end
valid.data=valid.data & isfinite(ybar) & isfinite(dhdy);
valid.combined=false(size(valid.data));
if ~any(valid.data); valid.combined=false(size(valid.data)); return; end

% omit segs with y vals too far from the median and recalculate.
valid.dist= abs(ybar-median(ybar(valid.data))) < params.ref_seg_Lsearch ;
Ymed=median(ybar(valid.dist & valid.data));
 
% select by y slope
% fit cross-track slope as a linear function of across-track distance, 
% delete pairs that are far from the curve. Iterate twice.
good=valid.dist & valid.data;
 
for k=1:2
    N_good=sum(good);
    if N_good==0; return; end
    this_degree= struct('x', double(diff(range(xpair(good)))>0),'y',double(diff(range(ybar(good)))>0));
        
    G_slope_y=[ones(N_good,1), ...
        poly2_fit_mat((xpair(good)-params.x_RGT_ctr)/y_scale, (ybar(good)-Ymed)/y_scale,   this_degree)];
        
    m_slope_y=(G_slope_y\dhdy(good));
    r_slope_y(good)=dhdy(good)-G_slope_y*(G_slope_y\dhdy(good));
    % changed 7/13/16 from iqr(r(good))
    sigma_slope_y= iqr(dhdy(good))/2;
    slope_threshold_y=max(0.01, max(3*sigma_slope_y, 3*median(sigma_dhdy(good))));
    good=isfinite(r_slope_y) & abs(r_slope_y-median(r_slope_y(good)))<slope_threshold_y;
end
valid.slope_y=good;

if DOPLOT
    axes( params.axes(1,1));
    plot(ybar, dhdy,'kx'); hold on;
    yvals=(min(D3a.y_RGT):5:max(D3a.y_RGT))';
    Gtemp=[ones(length(yvals),1),poly_fit_mat((yvals-Ymed)/y_scale, 1)];
    plot(yvals, Gtemp*m_slope_y,'b', yvals, Gtemp*m_slope_y+slope_threshold_y,'b--', yvals, Gtemp*m_slope_y-slope_threshold_y,'b--')
    plot(ybar(good_slope_y), dhdy(good_slope_y),'ro')
    ylabel('across-track slope');
end

% select by x slope (fit is by seg here);
good_other_criteria=ismember([D3a.rep, D3a.x_RGT], uRep_x0(valid.dist & valid.data, :),'rows') & isfinite(D3a.dh_fit_dx);
G_slope_x=[ones(size(good_other_criteria)), poly2_fit_mat((D3a.x_RGT-params.x_RGT_ctr)/y_scale, (D3a.y_RGT-Ymed)/y_scale, struct('x', 1,'y', 1))];
this_good=true(size(good_other_criteria));
for k=1:2
    these=this_good&good_other_criteria;
    these_cols=[true, max(G_slope_x(these, 2:end))-min(G_slope_x(these, 2:end)) > 0];
    m_slope_x=G_slope_x(these,these_cols)\D3a.dh_fit_dx(these);
    r_slope_x=D3a.dh_fit_dx-G_slope_x(these_cols)*m_slope_x;
    %sigma_slope_x=iqr(r_slope_x(this_good&good_other_criteria))/2;
    sigma_slope_x=iqr(D3a.dh_fit_dx(these))/2;
    slope_threshold_x=max([0.01, 3*sigma_slope_x, 3*median(D3a.sigma_dh_fit_dx(good))]);
    this_good=isfinite(r_slope_x) & abs(r_slope_x)< slope_threshold_x;
end
good_slope_x_seg=this_good;
% map the good (by seg) back to the pairs
% figure out which pairs have two good_slope_x values
valid.slope_x=false(size(uRep_x0,1),1);
seg_rep_x0=[D3a.rep, D3a.x_RGT];
for k=1:size(uRep_x0,1)
    these=seg_rep_x0(:,1)==uRep_x0(k,1) & seg_rep_x0(:,2)==uRep_x0(k,2);
    valid.slope_x(k)=sum(good_slope_x_seg(these))==2;
end

valid.Rep_x0=uRep_x0;
valid.combined=all([valid.dist, valid.slope_x, valid.slope_y, valid.data], 2);

if DOPLOT
    axes(params.axes(2,1));
    plot(D3a.y_RGT(good_dist_seg), D3a.dh_fit_dx(good_dist_seg),'kx'); hold on;
    plot(yvals, Gtemp*m_slope_x,'b', yvals, Gtemp*m_slope_x+slope_threshold_x,'b--', yvals, Gtemp*m_slope_x-slope_threshold_x,'b--')
    good=ismember(D3a.rep, good_reps);
    plot(D3a.y_RGT(good), D3a.dh_fit_dx(good),'ro');
    ylabel('along-track slope');
end
 
%-------------------------------------------------

function [P_sorted, deg_ind]=sort_poly2_components(P, deg, max_deg, deg_ind)
% return the polynomial coefficients, sorted first by x degree, then by x+y
% degree

if ~exist('deg_ind','var')
    [deg_x, deg_y]=meshgrid(0:max_deg.x, 0:max_deg.y);
    good_components=deg_x+deg_y > 0 & deg_x+deg_y <= max(max_deg.x, max_deg.y);
    deg_x=deg_x(:)';
    deg_x=deg_x(good_components);
    deg_y=deg_y(:)';
    deg_y=deg_y(good_components);
    [~, ind]=sort(deg_x+(deg_x+deg_y)/max(deg_x+deg_y+1));
    deg_x=deg_x(ind); deg_y=deg_y(ind);
    deg_ind=[deg_x(:), deg_y(:)];
end
% return if we're just looking for the polynomial coefficient order
if isempty(P); P_sorted=[]; return;end
    
[~, ind1, ind2]=intersect([deg.x(:), deg.y(:)], deg_ind,'rows');
P_sorted=zeros(1, size(deg_ind,1));
P_sorted(ind2)=P(ind1);

%-------------------------------------------------
function [G, coeff_degree]=poly2_fit_mat(x, y, degree)

% returns a matrix G such that G*m gives the values for a polynomial of
% degree (degree) at locations x,y for coefficients m.  Contains no
% zero-degree term
G=[];
[deg_x, deg_y]=meshgrid(0:degree.x, 0:degree.y);
good_components=deg_x+deg_y > 0 & deg_x+deg_y <= max(degree.x, degree.y);
deg_x=deg_x(:)';
deg_x=deg_x(good_components);
deg_y=deg_y(:)';
deg_y=deg_y(good_components);
Gx=real(exp(log(x(:))*deg_x));
Gy=real(exp(log(y(:))*deg_y));
Gx(~isfinite(Gx))=0;  Gx(:, deg_x==0)=1;
Gy(~isfinite(Gy))=0;  Gy(:, deg_y==0)=1;
G=Gx.*Gy;
coeff_degree.x=deg_x;
coeff_degree.y=deg_y;

%--------------------------------------------------
function varargout=reduce_by(fh, D, out_fields, index_field_name, field_vals)

index_field=D.(index_field_name);

for k=1:length(out_fields)
    varargout{k}=NaN(size(field_vals));
end

for k=1:length(field_vals)
    these=index_field==field_vals(k) & isfinite(D.x_RGT) & isfinite(index_field);
    if any(these)
        for kf=1:length(out_fields)
            if ~isfield(D,out_fields{kf}); 
                disp(['need ', out_fields{kf}]); 
                varargout{kf}(k)=NaN;
            else
                temp=D.(out_fields{kf})(these);
                if any(isfinite(temp))
                    varargout{kf}(k)=fh(temp(isfinite(temp)));
                else
                    varargout{kf}(k)=NaN;
                end
            end
        end
    end
end
 
%---------------------------------------------------
function varargout=regress_to(D, out_fields, index_field_name, index_field_ctr)

[varargout{1:nargout}]=deal(NaN);
G=[ones(numel(D.(index_field_name{1})),length(index_field_name)+1 )];
for k=1:length(index_field_name)
    G(:,k+1)= D.(index_field_name{k}) -index_field_ctr(k);
end
for k=1:length(out_fields)
    m=G\D.(out_fields{k});
    varargout{k}=m(1);
end

%--------------------------------------------------
function varargout=weighted_mean_by(D, weights, out_fields, index_field_name, field_vals)

index_field=D.(index_field_name);

for k=1:length(out_fields)
    varargout{k}=NaN(size(field_vals));
end

for k=1:length(field_vals)
    these=index_field==field_vals(k) & isfinite(D.x_RGT) & isfinite(index_field);
    if any(these) && sum(weights(these))>0
        
        for kf=1:length(out_fields)
            z=D.(out_fields{kf})(these);
            W=weights(these);
            ii=isfinite(z);
            if any(ii)
                varargout{kf}(k)=sum(W(ii).*z(ii))./sum(W(ii));
            end
        end
    end
    
end

%--------------------------------------------------
function [stats, summary]=assign_pass_quality_stats(N_reps, D3a, valid, out_fields)

% assign the output fields
for kf=1:length(out_fields)
    stats.(out_fields{kf})=NaN(1,N_reps);
end

rep_x0=[D3a.rep, D3a.x_RGT];
uRep=unique(rep_x0(:,1));

 % for the fit reps, find the best of the values
for k=1:length(uRep)
    this_rep=uRep(k);
    these=D3a.rep==uRep(k) & ismember(rep_x0, valid.Rep_x0(valid.combined,:),'rows');
    % if there are no pairs valid for the polynomial fit, find the seg that was reported  
    if ~any(these)
        these=D3a.rep==uRep(k) & valid.selected_segs;
    end
    if any(these)
        stats.ATL06_summary_zero_count(this_rep)=sum(D3a.ATL06_quality_summary(these)==0);
        stats.min_SNR_significance(this_rep)=min(D3a.SNR_significance(these));
        %stats.max_uncorr_reflectance=max(D3a.reflect_uncorr(these));
        stats.min_signal_selection_source(this_rep)=min(D3a.signal_selection_source(these));
    end
end

stats.mean_uncorr_reflectance=reduce_by(@mean, D3a, {'reflct_uncorr'},'rep', 1:N_reps);

summary=~(stats.min_signal_selection_source<=1 | stats.min_SNR_significance <0.02 | stats.ATL06_summary_zero_count >0);

%--------------------------------------------------
function PS=assign_pass_stats(D3a, weights, reps, fields)
 
for kf=1:length(fields)
    PS.(fields{kf})=NaN(size(reps(:)'));
end

%reduce_by(fh, D, out_fields, index_field_name, field_vals)
[PS.bsl_h_mean, PS.sigma_g_h, PS.sigma_g_x, PS.sigma_g_y, PS.laser_beam, PS.y_RGT,  PS.mean_pass_lat, PS.mean_pass_lon]=...
   reduce_by(@mean, D3a, {'bsl_h', 'sigma_g_h', 'sigma_g_x','sigma_g_y', 'laser_beam', 'y_RGT','lat_ctr','lon_ctr'}, 'rep', reps);
[PS.bslh_conf_best, PS.cloud_flg_best]=reduce_by(@min, D3a, {'bslh_conf', 'cloud_flg'},'rep', reps);
[ PS.h_robust_spread_mean, PS.reflct_uncorr_mean, PS.h_rms_mean, PS.tide_ocean_mean, PS.snr_mean] = ...
    weighted_mean_by( D3a, weights,  { 'h_robust_spread', 'reflct_uncorr', 'h_rms','tide_ocean', 'SNR'}, 'rep', reps);

 
  


    

