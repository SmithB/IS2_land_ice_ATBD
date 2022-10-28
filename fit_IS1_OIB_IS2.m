
params.E_RMS_d2z0dx2=2e-7;
params.E_RMS_d3zdx2dt=6e-8;

if ~exist('DEM','var')
    DEM=read_geotif('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_DEM/20140119_1453_10200100282DE200_102001002A716000-DEM_tr4x_LP100m.tif');
    
    DEM_img=DEM;
    DEM_img.x=DEM_img.x(1:5:end);
    DEM_img.y=DEM_img.y(1:5:end);
    DEM_img.z=DEM_img.z(1:5:end, 1:5:end);
    
    DEM_img.z=repmat(scale_to_byte(gradient(DEM.z), [-0.1 0.1]), [1 1 3]);

end
W=sqrt(diff(range(DEM.x).*diff(range(DEM.y))));
sigma=struct( 'smooth_dz',params.E_RMS_d3zdx2dt*W, 'smooth_z0', params.E_RMS_d2z0dx2*W, ...
    'line_spacing', 2.5e4);
if ~exist('D','var');
    ATL11_files=glob('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13/ATL11//run_1/*.h5');
    for k=1:length(ATL11_files); D11(k)=read_ATL11_h5(ATL11_files{k}); end
    clear temp DNR
    for k=1:length(D11); temp(k)=D11(k).non_repeat_data; end
    f=fieldnames(temp);
    for k=1:length(f)
        DNR.(f{k})=cat(1, temp.(f{k}));
    end
    DNR.x=ll2ps(DNR.lat, DNR.lon);
    
    clear temp D0 D
    for k=1:length(D11)
        temp(k)=D11(k).corrected_h;
    end
    f=fieldnames(temp);
    for k=1:length(f)
        D0.(f{k})=cat(1, temp.(f{k}));
    end
    good=D0.pass_quality_summary==1;
    D0.x=repmat(ll2ps(D0.ref_pt_lat, D0.ref_pt_lon), [1, size(D0.mean_pass_time, 2)]);
    D0=rmfield(rmfield(D0,'ref_pt_lat'),'ref_pt_lon');
    f=fieldnames(D0);
    for k=1:length(f)
        D1.(f{k})=D0.(f{k})(:);
        D1.(f{k})=D0.(f{k})(good(:));
    end
    
    clear D
    D.time=[datenum('September 12 2018')+D1.mean_pass_time; DNR.time];
    D.x=[D1.x; DNR.x];
    D.y=imag(D.x); D.x=real(D.x);
    D.z=[D1.pass_h_shapecorr; DNR.h_shapecorr];
    D.sigma=[D1.pass_h_shapecorr_sigma; DNR.h_shapecorr_sigma];
    D.sensor=[zeros(size(D1.mean_pass_time))+5; DNR.sensor];
    D=index_struct(D, isfinite(D.z));
end
time_nodes=[ datenum('Sep 15 2009'), datenum('Sep 15 2018')];

M_dz.XR=round_to(range(DEM.x), 1e3);
M_dz.YR=round_to(range(DEM.y), 1e3);
M_dz.dx=2e3;M_dz.dy=2e3;
grids_dz=fd_smooth_interp('define_grids', M_dz);

M_z0=M_dz; M_z0.dx=250;M_z0.dy=500;
grids_z0=fd_smooth_interp('define_grids', M_z0);

eqn_key=fd_smooth_interp('define_equation_types');

% z0 fit:
[G_z0, TOC_z0, ~,  D]=fd_smooth_interp('build_G_fit', D, grids_z0);
[Gc_z0, TOC_z0]=fd_smooth_interp('build_Gc', M_z0, grids_z0, TOC_z0);
% eliminate the flat rows
good=TOC_z0.eqn_type==eqn_key.smooth_z;
Gc_z0=Gc_z0(good,:);
TOC_z0.eqn_type=TOC_z0.eqn_type(good);
TOC_z0.eqn_id=TOC_z0.eqn_id(good);
TOC_z0.node_num=TOC_z0.node_num(good);


%space and time fit
[G_space, TOC_dz, ~,  D]=fd_smooth_interp('build_G_fit', D, grids_dz);
[Gc_space, TOC_dz]=fd_smooth_interp('build_Gc', M_dz, grids_dz, TOC_dz);
% eliminate the flat rows
good=TOC_dz.eqn_type==eqn_key.smooth_z;
Gc_space=Gc_space(good,:);
TOC_dz.eqn_type=TOC_dz.eqn_type(good);
TOC_dz.eqn_id=TOC_dz.eqn_id(good);
TOC_dz.node_num=TOC_dz.node_num(good);
 

% build the full fitting matrix:
%[Gdz     Gz0]
%[Gcdz    0
%[        Gcz0]

clear G_time;
%time node 1 and 2. t=zD - dzdt2(t2-t1) - dzdt1(t1-t)
els=D.time < time_nodes(1)  ;
G_time(els,1)=-(time_nodes(1)-D.time(els))/365.25;
G_time(els,2)=-(time_nodes(2)-time_nodes(1))/365.25;

els=D.time >=time_nodes(1) & D.time <= time_nodes(2);
G_time(els, 2)=-(time_nodes(2)-D.time(els))/365.25;

els=D.time > time_nodes(2);
G_time(els,3)=(D.time(els)-time_nodes(2))/365.25;

[GSc0.r, GSc0.c, GSc0.v]=find(Gc_space);
clear Gs Gsc TOC
N_x=size(G_space,2);
for kt=1:size(G_time,2)
    G1=G_space.*G_time(:, kt);
    [rr,cc,vv]=find(G1);
    GS(kt)=struct('r', rr, 'c', cc+(kt-1)*N_x, 'v', vv); 
    GSc(kt)=struct('r',(kt-1)*max(GSc0.r)+GSc0.r, 'c',(kt-1)*max(GSc0.c)+GSc0.c, 'v',GSc0.v);  
    TOC(kt).cols_dz=TOC_dz.cols.z+(k-1)*max(TOC_dz.cols.z);
    TOC(kt).dz_node=ones(size(TOC_dz.cols.z))*kt;
    TOC(kt).eqn_type=TOC_dz.eqn_type;
end

G_dz=sparse(cat(1, GS.r), cat(1, GS.c), cat(1, GS.v), size(G_time,1), size(G_space,2).*size(G_time,2));
Gc_dz=sparse(cat(1, GSc.r), cat(1, GSc.c), cat(1, GSc.v), max(cat(1, GSc.r)), size(G_space,2).*size(G_time,2));

G0=[[G_dz G_z0]; ...
    [Gc_dz, sparse([], [], [], size(Gc_dz,1), size(G_z0,2))]; ...
    [sparse([], [], [], size(Gc_z0,1), size(G_dz,2)), Gc_z0]];

TOC_all.eqn_type=[ones(size(G_z0,1),1); cat(1, TOC.eqn_type)+10; TOC_z0.eqn_type+20];

Cvals=zeros(size(TOC_all.eqn_type));
Cvals(TOC_all.eqn_type==1)=D.sigma.^2;
Cvals(TOC_all.eqn_type==10+eqn_key.smooth_z)=sigma.smooth_dz.^2/100000000;
Cvals(TOC_all.eqn_type==20+eqn_key.smooth_z)=sigma.smooth_z0.^2/100000;
  
mask=reshape(double(sum(G_space)>0), grids_dz.z.dims);
mask(mask==0)=NaN;
 

rows=true(size(G0,1),1);


for k=1:5
    m=my_lscov(G0(rows,:), [D.z(rows(1:length(D.z))); zeros(size(G0,1)-length(D.z),1)], 1./sqrt(Cvals(rows)));
    r=G0(1:length(D.z),:)*m-D.z;
    rs=r./D.sigma;
    
    sigma_hat=iqr(rs)/2;
    rows(1:length(D.z))=abs(r)<10 & abs(rs) < 3*sigma_hat;
    
    figure(1); clf;
    imagesc(mask.*reshape(m((1:prod(grids_dz.z.dims))+0*prod(grids_dz.z.dims)), grids_dz.z.dims))
    caxis([-10 10]);
    drawnow
end

clf
h=cheek_by_jowl(3,1, [0.1 0.1 0.8 0.8]);
axes(h(1));
image_struct(DEM_img); hold on;
imagesc(grids_dz.z.ctrs{2}, grids_dz.z.ctrs{1}, reshape(m((1:prod(grids_dz.z.dims))+0*prod(grids_dz.z.dims)), grids_dz.z.dims),'alphadata', 0.7*double(isfinite(mask)))
colormap(flipud(jet)); caxis([-10 10]);
these=D.time < time_nodes(1);
plot(D.x(these), D.y(these),'k.','markersize', 12);
these=these & rows(1:length(D.x));
plot(D.x(these), D.y(these),'w.','markersize', 6);

axes(h(2));
image_struct(DEM_img); hold on;
imagesc(grids_dz.z.ctrs{2}, grids_dz.z.ctrs{1}, reshape(m((1:prod(grids_dz.z.dims))+1*prod(grids_dz.z.dims)), grids_dz.z.dims),'alphadata', 0.7*double(isfinite(mask)))
colormap(flipud(jet)); caxis([-10 10]);
these=D.time > time_nodes(1) & D.time < time_nodes(2);
plot(D.x(these), D.y(these),'k.','markersize', 12);
these=these & rows(1:length(D.x));
plot(D.x(these), D.y(these),'w.','markersize', 6);

axes(h(3));
image_struct(DEM_img); hold on;
imagesc(grids_dz.z.ctrs{2}, grids_dz.z.ctrs{1}, reshape(m((1:prod(grids_dz.z.dims))+2*prod(grids_dz.z.dims)), grids_dz.z.dims),'alphadata', 0.7*double(isfinite(mask)))
colormap(flipud(jet)); caxis([-10 10]);
these=D.time > time_nodes(2);
plot(D.x(these), D.y(these),'k.','markersize', 12);
these=these & rows(1:length(D.x));
plot(D.x(these), D.y(these),'w.','markersize', 6);













