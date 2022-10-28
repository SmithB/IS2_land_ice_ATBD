

% plots showing ATL14/15 output accuracy

%[~, out]=unix('ls v9_hdf/PIG_ATL11_dh_clouds/run*/ATL13_solution.mat');
[~, out]=unix('ls /Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13/ATL14/run_*/ATL13_solution.mat');
% [~, out]=unix('ls v9_hdf/PIG_ATL11_with_dh_run_*/tau=2.0//ATL13_solution.mat');
out=strsplit(deblank(out));
for k=1:length(out)
    temp=load(out{k},'dzbar');
    temp.dzbar.sigma=sqrt(mean((temp.dzbar.z_est(:,:,1:12)-temp.dzbar.z_true(:,:,1:12)).^2,3));
    dzbar(k)=temp.dzbar;
end;
sig=cat(3, dzbar.sigma);

sigma_all=(sqrt(mean(sig.^2,3)));

filename='Dh_real_and_est.gif';
delete(filename);
clf; set(gcf,'position', [550 35 1090 1000],'color', [1 1 1],'defaultaxesfontsize', 14,'defaultaxesfontname','times');
h=cheek_by_jowl(3,1, [0.05 0.07 0.8 0.9]);
dzg=PIG_dh_function(xg, yg, tg*365);
dzg=dzg-repmat(dzg(:,:,1), [1 1 size(dzg,3)]);
dz.z=dz.z-repmat(dz.z(:,:,1), [1 1 size(dz.z,3)]);

MM=double(dz.mask); MM(dz.mask==0)=NaN;
hb=axes('position', [0.86 0.39 0.02 0.58]);
axes(hb); imagesc([0 1], [-11.5:.1:11.5], [-11.5:.1:11.5]'); caxis([-11.5 11.5]);set(hb,'yaxislocation','right','xtick', []); ylabel(hb,'dh WRT h(0), m');
set(gca,'ylim', [-10.5 10.5]);
colormap([0 0 0; flipud(dzdt_cpt(128))]);
set(h(1:2),'visible','off');
rc1=[14, 42];
rc2=[14, 27];
for k=1:13;
    axes(h(1)); cla
    imagesc(dz.x, dz.y, (dz.z(:,:,k)).*MM); caxis([-11.5 11.5]); axis xy equal tight;
    XL=get(gca,'xlim'); YL=get(gca,'ylim');
    ht=text(XL*[0.95 0.05]', YL*[0.1 0.9]', sprintf('T=%d months', round(dz.t(k)*365/30.4)),'color', [1 1 1],'fontweight','bold','fontsize', 14);
    hold on;plot(dz.x(rc1(2)), dz.y(rc1(1)),'wo','markersize', 8,'linewidth', 2);
    plot(dz.x(rc2(2)), dz.y(rc2(1)),'wo','markersize', 8, 'linewidth', 2);
    els=abs(D.t-dz.t(k)*365)<.25*365/2;
    
    plot(unique(D.x(els)+1i*D.y(els)),'k.','markersize', 1)
    axes(h(2)); cla;
    imagesc(dz.x, dz.y, (dzg(:,:,k)-dzg(:,:,1)).*MM); caxis([-11.5 11.5]); axis xy equal tight;
    hold on;plot(dz.x(rc1(2)), dz.y(rc1(1)),'wo','markersize', 8,'linewidth', 2);
    plot(dz.x(rc2(2)), dz.y(rc2(1)),'wo','markersize', 8, 'linewidth', 2);
    plot(unique(D.x(els)+1i*D.y(els)),'k.','markersize', 1)
    set(h(1:2),'visible','off');
    set(findobj(gca,'type','line','marker','o'),'markerfacecolor', [0 0 0]);
    
    
    axes(h(3)); cla; hold on
    plot( (dz.t*365/30.4), squeeze(dzg(rc1(1), rc1(2), :)),'ko');
    plot( (dz.t(1:k)*365/30.4), squeeze(dz.z(rc1(1), rc1(2), 1:k)),'b','linewidth', 3);
    plot( (dz.t*365/30.4), squeeze(dzg(rc2(1), rc2(2), :)),'ko');
    plot( (dz.t(1:k)*365/30.4), squeeze(dz.z(rc2(1), rc2(2), 1:k)),'r','linewidth', 3);
    ylabel('dh WRT h(0), m');xlabel('month');
    set(gca,'xlim', [-0.5 36.5] ,'ylim', [-10.5 10.5]);
    grid on;
    
    drawnow; frame=getframe(gcf);
    im=frame2im(frame);
    [A, map]=rgb2ind(im, 256);
    if k==1;
        imwrite(A, map, filename,'gif', 'LoopCount', Inf,'DelayTime', 1);
    else
        imwrite(A, map, filename,'gif', 'DelayTime', 1,'WriteMode', 'append');
    end
end

imwrite(A, map, filename,'gif', 'DelayTime', 5,'WriteMode', 'append');



% spin the DEMs movie;
II=read_geotif('PIG_DEM/20140119_1453_10200100282DE200_102001002A716000-DEM_tr16x.tif');
II.z(II.z==0)=NaN;
M0=double(z0.mask); M0(z0.mask==0)=NaN;
ha=axes('position', [0 0 1 1]);
hold on;
hs(1)=surf(II.x(1:10:end), II.y(1:10:end), double(II.z(1:10:end, 1:10:end)));shading interp; material dull;  hs(2)==surf(z0.x, z0.y-3e4, full(z0.z.*M0)); shading interp; material dull; lighting phong;  hl=light;
set(gca,'dataaspectratio', [ 1 1 .01]);
set(gca,'visible','off')
view(0, 70);
filename='two_DEM_movie.gif';
for theta=0:2:360;
    view( theta, 70);
    drawnow; frame=getframe(gcf);
    im=frame2im(frame);
    [A, map]=rgb2ind(im, 256);
    if theta==0;
        imwrite(A, map, filename,'gif', 'LoopCount', Inf,'DelayTime', 1/8);
    else
        imwrite(A, map, filename,'gif', 'DelayTime', 1/8,'WriteMode', 'append');
    end
end

imwrite(A, map, filename,'gif', 'DelayTime', 5,'WriteMode', 'append');


% make an error map
[II.gx]=gradient(II.z, II.x, II.y);
figure(123);clf
imagesc(II.x(1:5:end), II.y(1:5:end), repmat(scale_to_byte(II.gx(1:5:end, 1:5:end), [-.05 0.05]), [1 1 3]),'alphadata', double(isfinite(II.gx(1:5:end, 1:5:end))));
axis xy equal tight
hold on;
plot(unique(D.x+1i*D.y),'k.','markersize', 1);
imagesc(dzbar(1).x, dzbar(1).y, sigma_all,'alphadata', 0.7*isfinite(sigma_all));
hb=colorbar;
ylabel(hb,'RMSE(\deltaz), m');
set(gca,'visible','off');
figure(124)
hold on;
for k=1:10; plot(dz.t, squeeze(dzbar(k).z_est(2, 3,:))); end
plot(dz.t, squeeze(dzbar(1).z_true(2,3,:)),'ko');
xlabel('year'); ylabel('\delta h, m');


 