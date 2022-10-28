
ATL11_files=(glob('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13/ATL11//run_1/*.h5'));

clear D6;
for k=1:length(ATL11_files)
    [~, D6(k)]=read_ATL11_h5(ATL11_files{k});
end

xx=cat(1, D6.x_PS_ctr); 
yy=cat(1, D6.y_PS_ctr);
tt=cat(1, D6.time);

[xg,yg]=meshgrid(DEM.x(1:50:end), DEM.y(1:50:end));
aa=double(isfinite( DEM.z(1:50:end, 1:50:end)));

clf;
hax=cheek_by_jowl(2,1,[0.1 0.2 0.8 0.6]);

axes(hax(1)); 
imagesc(DEM.x, DEM.y, repmat(scale_to_byte(DEM.gx, [-.1 .1]), [1 1 3])); axis xy equal tight
hold on;
hi=imagesc(xg(1,:), yg(:,1), PIG_dh_function(xg, yg, zeros(size(xg))),'alphadata', aa*.125);

alpha_map.alpha=[1 1 .25 .25 1 1];
alpha_map.z=    [-100 -2 -.25 .25 2 100];

axes(hax(2));
imagesc(DEM.x, DEM.y, repmat(scale_to_byte(DEM.gx, [-.1 .1]), [1 1 3])); axis xy equal tight
hold on;

axes('position', [.2 .1 .6 .025]);
imagesc([-10 10], [0 1], [-10 :10]);
set(gca,'ytick', [],'fontsize', 13); 
xlabel('\delta z, m');

set(hax,'xtick', [],'ytick', []);
t_vals=0:15:(3*360);
clear F
for k=1:length(t_vals);
    this_dz=PIG_dh_function(xg, yg, t_vals(k)*ones(size(xg)));
    this_alpha=interp1(alpha_map.z, alpha_map.alpha, this_dz);
    this_alpha(~isfinite(this_alpha))=0;
    axes(hax(1));
    caxis([-10 10]);
    set(hi,'cdata',this_dz,'alphadata', aa.*this_alpha*.5);
    title(hax(1), sprintf('%d days after launch', t_vals(k)));
    
    axes(hax(2));
    delete(findobj(gca,'type','line'));
    els=abs(tt-t_vals(k))<7.5;
    plot_colored_points(xx(els)+1i*yy(els), PIG_dh_function(xx(els), yy(els), tt(els)), [-10:.25:10], [], get(gcf,'colormap'), true);
    set(hax(2),'xlim', range(DEM.x),'ylim', range(DEM.y));
    drawnow;
    F(k)=getframe(gcf);
end
V=VideoWriter('sim_dh_map.avi');
V.FrameRate=4;
V.open
for k=1:length(F);
    writeVideo(V,F(k));
end
close(V);


