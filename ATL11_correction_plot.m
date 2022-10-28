%------------------------------------------------------------------
function ATL11_correction_plot(D11, D6, seg,  params)

x0_ATC=D11.reference_point.x_ATC(seg);
y0_ATC=D11.reference_point.y_ATC(seg);

temp=index_struct(D6, abs(D6.x_RGT-x0_ATC)<=60);
temp.dx_RGT=temp.x_RGT-x0_ATC;
temp.dy_RGT=temp.y_RGT-y0_ATC;

[G0, coeff_degree]=ATL11_proc_ATBD_Nseg('poly2_fit_mat',temp.dx_RGT(:)/100, temp.dy_RGT(:)/100, struct('x', 4,'y', 3));
G=zeros(numel(temp.dx_RGT), size(params.poly_exp_x,2));
for kC=1:size(params.poly_exp_x',2)
    this_col=coeff_degree.x==params.poly_exp_x(kC) & coeff_degree.y==params.poly_exp_y(kC);
    G(:, kC)=G0(:,this_col);
end
temp.poly_corr=G*(D11.ref_surf.poly_ref_surf(seg, :)');
C_m=diag(D11.ref_surf.poly_ref_surf_sigma(seg,:).^2);
temp.poly_corr_sigma=sqrt(diag(G*(C_m*G')));

sigma_hc=sqrt(temp.h_LI_sigma.^2+temp.poly_corr_sigma.^2);
% sigma_hc_plus=sqrt(temp.sigma.^2+temp.poly_corr_sigma.^2+temp.bias_50m.^2);
% if all(isfinite(sigma_hc_plus))
%     best=find(sigma_hc_plus==min(sigma_hc_plus));
% else
%     best=find(sigma_hc==min(sigma_hc));
% end
% temp.h_shapecorr_sigma=sigma_hc;
% temp.h_shapecorr=temp.z-temp.poly_corr;
 

% now make the correction surface
[dx, dy]=meshgrid([-120:10:120], [-200:10:200]);

[G0, coeff_degree]=ATL11_proc_ATBD_Nseg('poly2_fit_mat',dx(:)/100, dy(:)/100, struct('x', 4,'y', 3));
G=zeros(numel(dx), size(params.poly_exp_x,2));
for kC=1:size(params.poly_exp_x',2)
    this_col=coeff_degree.x==params.poly_exp_x(kC) & coeff_degree.y==params.poly_exp_y(kC);
    G(:, kC)=G0(:,this_col);
end
poly_corr_surf=reshape(G*(D11.ref_surf.poly_ref_surf(seg, :)'), size(dx));
C_m=diag(D11.ref_surf.poly_ref_surf_sigma(seg,:).^2);
poly_corr_sigma=reshape(sqrt(diag(G*(C_m*G'))), size(dx));

z0=median(D11.corrected_h.pass_h_shapecorr(seg, D11.corrected_h.pass_included_in_fit(seg,:)==1));

hax(1)=subplot(2,1,1); cla
set(gca,'position', get(gca,'position').*[1 1 .9 1]);
[h, h_bar(1)]=plot3_colored_points(temp.x_RGT-x0_ATC+1i*(temp.y_RGT-y0_ATC), temp.h_LI, temp.time/365, [0:.1:3]);
hold on;
% mark the low-quality points
non_sel_pts=~ismember(temp.rep, find(D11.corrected_h.pass_included_in_fit(seg,:)==1));
plot3(temp.x_RGT(non_sel_pts)-x0_ATC, temp.y_RGT(non_sel_pts)-y0_ATC, temp.h_LI(non_sel_pts),'o','color', [0.6 0.6 0.6],'linewidth', 2);
set(gca,'dataaspectratio', [1 1 .03]);
ylabel(h_bar(1),'time after launch, years');
set(h(ishandle(h)),'markersize', 16);
hold on;
hs=surf(dx, dy, z0+poly_corr_surf); set(hs, 'facecolor', [0.6 0.6 0.6]);
hs1=surf(dx, dy, z0+poly_corr_surf+poly_corr_sigma);set(hs1,'facecolor','none','edgecolor','k');
hs2=surf(dx, dy, z0+poly_corr_surf-poly_corr_sigma); set(hs2, 'facecolor','none','edgecolor','k');
xlabel('along-track, m'); 
ylabel('across-track, m');
zlabel('height, m');
grid on
view([-56 10]);


hax(2)=subplot(2,1,2);
set(gca,'position', get(gca,'position').*[1 1 .9 1]);

iii=sub2ind(size(D11.corrected_h.pass_h_shapecorr), seg*ones(size(temp.rep)), temp.rep);
[h, h_bar(2)]=plot3_colored_points(temp.x_RGT-x0_ATC+1i*temp.y_RGT-y0_ATC, temp.h_LI-D11.corrected_h.pass_h_shapecorr(iii), temp.time/365, [0:.1:3]);
set(h(ishandle(h)),'markersize', 16);

set(gca,'zlim', diff(get(hax(1),'zlim'))*[-.5 .5]);

ylabel(h_bar(2),'time after launch, years');

xlabel('along-track, m'); 
ylabel('across-track, m');
zlabel('height correction, m');
set(gca,'dataaspectratio', [1 1 .03])
grid on
view([-56 10]);
set(hax, 'xlim', range(dx(1,:)),'ylim', range(dy(:,1)));

figure; 
errorbar(D11.corrected_h.mean_pass_time(seg,:)/365, D11.corrected_h.pass_h_shapecorr(seg,:), sqrt(D11.corrected_h.pass_h_shapecorr_sigma(seg,:).^2+D11.corrected_h.pass_h_shapecorr_sigma_systematic(seg,:).^2),'ko')
hold on;
NotFit=D11.corrected_h.pass_included_in_fit(seg,:)==0;
plot(D11.corrected_h.mean_pass_time(seg,NotFit)/365, D11.corrected_h.pass_h_shapecorr(seg,NotFit),'rx','linewidth', 3);

xlabel('time after launch, yr'); ylabel('corrected h, m');

% make a movie
if false
clear f;
dtheta=0:2:720;
set([hax(2); h_bar(2); get(h_bar(2),'children'); findobj(hax(1),'type','surface'); findobj(hax(2),'type','line')],'visible','off');
for k=1:length(dtheta)/2
    axes(hax(1)); view([-56-dtheta(k), 10]);
    axes(hax(2)); view([-56-dtheta(k), 10]);
    f(k)=getframe(gcf);
end
set([hax(2); h_bar(2);  get(h_bar(2),'children'); findobj(hax(1),'type','surface'); findobj(hax(2),'type','line')],'visible','on');
for k=length(dtheta)/2+0.5:length(dtheta)
    axes(hax(1)); view([-56-dtheta(k), 10]);
    axes(hax(2)); view([-56-dtheta(k), 10]);
    f(k)=getframe(gcf);
end

! rm ThreeD_ATL11.avi
V=VideoWriter('ThreeD_ATL11.avi')
V.FrameRate=15;
open(V);
for k=1:length(f)
    writeVideo(V,f(k));
end

close(V);




end
 