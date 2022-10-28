 function varargout=add_OIB_to_ATL11(varargin)

if nargout>0
    [varargout{1:nargout}]=feval( varargin{:});
else
    feval(varargin{:});
end

%------------------------------------------------------------------
function [D11, D200]=generate_corrections(filename)

disp(filename);
[D11, ~, params_11]=read_ATL11_h5(filename);
D11.corrected_h.mean_pass_time=D11.corrected_h.mean_pass_time+datenum('september 12 2018');
[D200, rgt_xhat, rgt_yhat, params_OIB]=add_OIB_to_ATL11('read_OIB_data', D11);
N_pts=size(D11.reference_point.rgt_azimuth,1);
clear D_corr temp;
for k=1:N_pts
    temp=add_OIB_to_ATL11('correct_one_point', D11, D200(k), k, rgt_xhat(k,:), rgt_yhat(k,:), params_11);
    if isempty(temp); continue; end
    for kk=1:length(temp)
        temp(kk).ref_pt=zeros(size(temp(kk).time))+D11.reference_point.seg_count(k);
    end
    D_corr(k)=temp;
    
    % and plot....
    if false %~isempty(D_corr(k).z)
        figure(10+k);
        clf
        good=D11.corrected_h.pass_quality_summary(k,:)==0;
        errorbar(D11.corrected_h.mean_pass_time(k,good), D11.corrected_h.pass_h_shapecorr(k,good), sqrt(D11.corrected_h.pass_h_shapecorr_sigma(k,good).^2+D11.corrected_h.pass_h_shapecorr_sigma_systematic(k,good).^2),'o','color', [0 0.7 0]);
        hold on;
        errorbar(D_corr(k).time, D_corr(k).h_shapecorr, D_corr(k).h_shapecorr_sigma,'x','color', [0.7 0 0.2]);
    end
    
    
end

f=fieldnames(D_corr); 

f=f(~ismember(f, {'EN', 'x','y','slope_x','slope_y','noise_50m','bias_50m'}));
for kf=1:length(f)
    D11.non_repeat_data.(f{kf})=cat(1, D_corr.(f{kf}));
end

 
%------------------------------------------------------------------
function D_corr=correct_one_point(D11, D200, seg, rgt_xhat, rgt_yhat, params)

D200.EN=platte_carre_WGS84([D200.lat, D200.lon],  [D11.corrected_h.ref_pt_lat(seg) D11.corrected_h.ref_pt_lon(seg)],'forward');
D200.dx_RGT=D200.EN(:,1)*rgt_xhat(1)+D200.EN(:,2)*rgt_xhat(2);
D200.dy_RGT=D200.EN(:,1)*rgt_yhat(1)+D200.EN(:,2)*rgt_yhat(2);


day=floor(D200.time);
days=unique(day);
for kD=1:length(days)
    temp=index_struct(D200, day==days(kD));
    [G0, coeff_degree]=ATL11_proc_ATBD_Nseg('poly2_fit_mat',temp.dx_RGT(:)/100, temp.dy_RGT(:)/100, struct('x', 4,'y', 3,'yscale', 100));
    G=zeros(numel(temp.dx_RGT), size(params.poly_exp_x,2));
    for kC=1:size(params.poly_exp_x',2)
        this_col=coeff_degree.x==params.poly_exp_x(kC) & coeff_degree.y==params.poly_exp_y(kC);
        G(:, kC)=G0(:,this_col);
    end
    temp.poly_corr=G*(D11.ref_surf.poly_ref_surf(seg, :)');
    C_m=diag(D11.ref_surf.poly_ref_surf_sigma(seg,:).^2);
    temp.poly_corr_sigma=sqrt(diag(G*(C_m*G')));
    
    % add xy errors for IS2, etc
    temp.sigma_geo_xy=zeros(size(temp.z));
    temp.sigma_geo_xy(temp.sensor < 1)=15;
    temp.sigma_geo_xy(temp.sensor >=1)=0.3;
    % temporary estimate of surface slope
    slope_mag=max(abs(temp.slope_x+1i*temp.slope_y), abs(D11.ref_surf.fit_X_slope(seg)+1i*D11.ref_surf.fit_Y_slope(seg)));
    
    sigma_hc=sqrt(temp.sigma.^2+temp.poly_corr_sigma.^2 + (slope_mag.*temp.sigma_geo_xy).^2);
    sigma_hc_plus=sqrt(temp.sigma.^2+temp.poly_corr_sigma.^2+temp.bias_50m.^2 + (slope_mag.*temp.sigma_geo_xy).^2);
    if all(isfinite(sigma_hc_plus)) 
        best=find(sigma_hc_plus==min(sigma_hc_plus));
    else
        best=find(sigma_hc==min(sigma_hc));
    end
    temp.h_shapecorr_sigma=sigma_hc;
    temp.h_shapecorr=temp.z-temp.poly_corr;
    temp.slope_mag=sqrt(temp.slope_x.^2+temp.slope_y.^2);
    D_corr(kD)=index_struct(temp, best);
    D_corr(kD).EN=temp.EN(best,:);
end
if ~exist('D_corr','var'); D_corr=[]; return; end

ff=fieldnames(D_corr); 
for kf=1:length(ff); 
    D_corr(1).(ff{kf})=cat(1, D_corr.(ff{kf}));
end
D_corr=D_corr(1);

%------------------------------------------------------------------------- 
function [D200, rgt_xhat, rgt_yhat, params]=read_OIB_data(D11, params)
if ~exist('params','var')
    params=get_master_indices;
end

xy=ll2ps(D11.corrected_h.ref_pt_lat, D11.corrected_h.ref_pt_lon);
delta_x=([-1 0 1 -1 0 1 -1 0 1]+1i*[-1 -1 -1 0 0 0 1 1 1])*250; 
bins=unique(round_to(xy(:)+delta_x, 1e4));

fields={'x','y','time','z','sensor', 'slope_x', 'slope_y','noise_50m','bias_50m'};
clear D_OIB;
D=index_point_data_h5('read_from_index', bins, params.master_index_IS, [fields, 'reflctUC', 'IceSVar', 'satElevCorr', 'gval_rcv']);

% run standard parameter-based cleanup
Thresholds=select_filter_params;
IS_els=find(D.sensor < 1);
good=false(size(IS_els));
period=round(D.sensor(IS_els)*100);
uP=unique(period);
Thresholds.reflctUncorr=max(Thresholds.reflctUncorr, 0.05);
for kP=1:length(uP);
    p_els=IS_els(period==uP(kP));
    good(period==uP(kP))=D.reflctUC(p_els) > Thresholds.reflctUncorr(uP(kP)) & ...
        D.IceSVar(p_els) < Thresholds.IceSVar(uP(kP)) & D.satElevCorr(p_els) < 0.25  & D.gval_rcv(p_els) >=13;
    if (uP(kP)) < 5
        good(period==uP(kP))=good(period==uP(kP)) & D.gval_rcv(p_els) < 250;
    end
end
delete_els=IS_els(~good);
good=true(size(D.sensor));
good(delete_els)=false;
D.z=D.z+D.satElevCorr;
D_OIB(1)=index_struct(D,good, fields);

D_OIB(2)=index_point_data_h5('read_from_index', bins, params.master_index_Qfit, fields);
D_OIB(3)=index_point_data_h5('read_from_index', bins, params.master_index_LVIS, fields);

clear temp
for kF=1:length(fields)
   temp.(fields{kF})=cat(1, D_OIB.(fields{kF})); 
end
D_OIB=temp;

[D_OIB.lat, D_OIB.lon]=ps2ll(D_OIB.x, D_OIB.y);
 
N_pts=size(D11.reference_point.rgt_azimuth,1);
rgt_xhat=[sind(D11.reference_point.rgt_azimuth) cosd(D11.reference_point.rgt_azimuth)];
rgt_yhat=[-cosd(D11.reference_point.rgt_azimuth) sind(D11.reference_point.rgt_azimuth)];

pt0=ceil(N_pts/2);

lat0=D11.corrected_h.ref_pt_lat(pt0);
lon0=D11.corrected_h.ref_pt_lon(pt0);
xhat_0=rgt_xhat(pt0,:);
yhat_0=rgt_yhat(pt0,:);
x_ATC_0=D11.reference_point.x_ATC(pt0);
 
EN_11=platte_carre_WGS84([D11.corrected_h.ref_pt_lat, D11.corrected_h.ref_pt_lon], [lat0, lon0], 'forward');

EN_OIB=platte_carre_WGS84([D_OIB.lat, D_OIB.lon], [lat0, lon0], 'forward');
D_OIB.EN=EN_OIB;
y_local=yhat_0(1)*EN_OIB(:,1)+yhat_0(2)*EN_OIB(:,2);
x_local=x_ATC_0+xhat_0(1)*EN_OIB(:,1)+xhat_0(2)*EN_OIB(:,2);
D_OIB_L=index_struct(D_OIB, abs(y_local)<500);
 
els=D_OIB_L.sensor>0;
D_OIB_L.sigma=0.1+zeros(size(D_OIB_L.time));
D_OIB_L.sigma(els)=sqrt(0.03^2+D_OIB_L.noise_50m(els).^2);
D_OIB_L.sigma(D_OIB_L.sensor<1)=0.15; 
% now loop over each reference point, collect the data within 200 m
for kP=1:N_pts
    els=(D_OIB_L.EN(:,1)-EN_11(kP,1)).^2+(D_OIB_L.EN(:,2)-EN_11(kP,2)).^2 < 200^2;
    D200(kP)=index_struct(D_OIB_L, els);
end

if false
    clf;
    plot(EN(:,1), EN(:,2),'ko'); hold on;
    plot((EN(:,1)+[0 50].*rgt_xhat(:,1))', (EN(:,2)+[0 50].*rgt_xhat(:,2))');
    plot((EN(:,1)+[0 50].*rgt_yhat(:,1))', (EN(:,2)+[0 50].*rgt_yhat(:,2))');
    
    axis xy equal
end

%---------------------------------------------------------------------------

function params=get_master_indices


params.h5_source.IS='/Volumes/insar7/gmap/oib_database/glas/AA/rel_634/master_index_h5.mat';
params.h5_source.Reigl='/Volumes/insar7/gmap/oib_database/riegl/AA/AA_srfelev.h5';
params.h5_source.ATM='/Volumes/insar7/gmap/oib_database/ATM_Qfit/Antarctica/master_index_h5.mat';
params.h5_source.LVIS='/Volumes/insar6/gmap/oib_database/LVIS/Antarctica/master_index.mat';

temp=load(params.h5_source.IS);
params.master_index_IS=temp.master_index;
params.Reigl_h5_file=params.h5_source.Reigl;
temp=load( params.h5_source.ATM);
params.master_index_Qfit=temp.master_index;
temp=load(params.h5_source.LVIS);
params.master_index_LVIS=temp.master_index;



%----------------------------------
function correct_files(filenames)

 
for k0=1:length(filenames)
    
    D11_file=filenames{k0};
    [data_dir, D11_base]=fileparts(filenames{k0});
     
    if ~exist([data_dir,'/backup'],'dir'); mkdir([data_dir,'/backup']); end
    if ~exist([data_dir,'/backup/', D11_base,'.h5'],'file')
        unix([ 'cp ',filenames{k0},' ', data_dir,'/backup']);
    end
    
    [D11, D200]=generate_corrections(filenames{k0});
    ff=fieldnames(D11.non_repeat_data);
    for kf=1:length(ff)
        this_field=['/non_repeat_data/', ff{kf}];
        temp=getfield(D11.non_repeat_data, ff{kf});
        if ~isempty(temp)
            try
                h5create(D11_file, this_field, [inf, size(temp,2)],'ChunkSize', [min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
            catch
            end
            h5write(D11_file, this_field,  double(temp), [1,1], size(temp));
        end
    end
        
    D200_dir=[data_dir,'/D200'];
    if ~exist(D200_dir,'dir'); mkdir(D200_dir);end

    D200_file=[D200_dir,'/', D11_base,'.h5'];
    
    ff=fieldnames(D200);
    for kf=1:length(ff)
        this_field=['/',ff{kf}];
        temp=cat(1, D200.(ff{kf}));
        try
            h5create(D200_file, this_field, [inf, size(temp,2)],'ChunkSize', [ min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
        catch
        end
        h5write(D200_file, this_field,  double(temp), [1 1], size(temp));
    end
end










