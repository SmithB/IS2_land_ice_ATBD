function [D_06, dh_hist, dhR]=atl06_dh_plot(thefile, D3)
fprintf(1,'working on %s', thefile);

D_06=read_ASAS_ATL06(thefile);

%this_cmap=[summer; flipud(autumn)]*.7;
%this_cmap=my_rgb_cpt(128)*.7+.1;
this_cmap=jet*.7+.1;

for kp=1:3;
    cloud_rat=D_06(kp).n_fit_photons./D_06(kp).w_surface_window_final;
    good_cloud=[cloud_rat(:,1)>1, cloud_rat(:,2) >4];
    good_AT=ATL06_AT_filter(D_06(kp), 2);
    qs=D_06(kp).h_robust_sprd>1 | D_06(kp).h_li_sigma > 1  | D_06(kp).snr_significance > 0.02;
    good=good_cloud & good_AT & (qs==0);
    D_06(kp).h_li(good==0)=NaN;
end
    

 
latR=range(D_06(2).latitude(isfinite(D_06(2).latitude)));
if max(latR) > 0
    hemisphere=1;
else
    hemisphere=-1;
end


MI_Q=getappdata(0,'master_index_Qfit');
if isempty(MI_Q);
    load /Volumes/insar7/gmap/oib_database/ATM_Qfit/Antarctica/master_index_h5.mat
    MI_Q=master_index;
    setappdata(0,'master_index_Qfit', MI_Q);
end
MI_L=getappdata(0,'master_index_LVIS');
if isempty(MI_L);
    load /Volumes/insar5/gmap/OIB_data/LVIS/AA/master_index_h5.mat
    MI_L=master_index;
    setappdata(0,'master_index_LVIS', MI_L);
end
 

for kP=1:3;
    D_06(kP).latitude(D_06(kP).latitude==0)=NaN;
end

for kP=1:3
    if hemisphere==1
        D_06(kP).x=gl_ll2ps(D_06(kP).latitude, D_06(kP).longitude);
    else
        D_06(kP).x=ll2ps(D_06(kP).latitude, D_06(kP).longitude);
    end
    
    browse_axes=findobj('tag','browse_axes');
    XL=get(browse_axes,'xlim');
    YL=get(browse_axes,'ylim');
    xx=D_06(kP).x;
    els=real(xx) > XL(1) & real(xx) < XL(2) & imag(xx) > YL(1) & imag(xx) < YL(2);
    els=any(els,2);
    D_06(kP)=index_struct(D_06(kP), els);
    D_06(kP).time=D_06(kP).delta_time/24/3600+datenum('jan 1 2018');
    D_06(kP).year=(D_06(kP).time-datenum('jan 1 2010'))/365.25+2010;
end
 
xx=cat(1, D_06.x);
XR=range(real(xx(:)));
YR=range(imag(xx(:)));
XC=mean(XR);
YC=mean(YR);
HW=max(diff(XR), diff(YR))/2;

if hemisphere==1
    MOS='/Volumes/ice1/ben/MOG/2005/mog100_2005_hp1_v1.1.tif';
else   
    MOS='/Volumes/ice1/ben/MOA/2012/moa125_hp1_2013310_2014067_v33.tif';
end
II=read_geotif_xy(MOS, XC+[-HW HW], YC+[-HW HW]*0.7);
figure(4); clf; set(gcf,'defaultaxesfontsize', 16,'inverthardcopy','off','color','w','units','inches'); 
pos=get(gcf,'position'); 
set(gcf,'papersize', pos(3:4)+0.5);

hax(1)=axes('position',[0.15 0.6 0.8 0.3]);
imagesc((II.x-XL(1))/1000, II.y/1000, II.z); axis xy equal
colormap(gray); hold on; 

for kP=1:3

        [xg, yg]=grid_box_cover(XL, YL, 1e4, true);
        D_ATM=index_point_data_h5('read_from_index', xg+1i*yg, MI_Q, {'x','y','z', 'time','sensor'});
        D_LVIS=index_point_data_h5('read_from_index', xg+1i*yg, MI_L, {'x','y','z', 'time','sensor'});
        
        D_ATM=index_struct(D_ATM,  D_ATM.x>II.x(1) & D_ATM.x < II.x(end) & D_ATM.y > II.y(1) & D_ATM.y < II.y(end)); 
        D_LVIS=index_struct(D_LVIS, D_LVIS.x>II.x(1) & D_LVIS.x < II.x(end) & D_LVIS.y > II.y(1) & D_LVIS.y < II.y(end)); 
                
        D_ATM.year=(D_ATM.time-datenum('jan 1 2010'))/365.25+2010;
        D_LVIS.year=(D_LVIS.time-datenum('jan 1 2010'))/365.25+2010;
end
 
axes(hax(1)); 
L=(D_ATM.x-XL(1))/1000;
plot_colored_points(L+1i*D_ATM.y/1000, D_ATM.year, 2002:2020, [], this_cmap, true);
L=(D_LVIS.x-XL(1))/1000;
plot_colored_points(L+1i*D_LVIS.y/1000, D_LVIS.year, 2002:2020, [], this_cmap, true);
L=(real(cat(1, D_06.x))-XL(1))/1000;
plot_colored_points(L+1i*imag(cat(1, D_06.x))/1000, cat(1, D_06.year), 2002:2020, [], this_cmap, true);

setappdata(gca,'D_ATM', D_ATM);
setappdata(gca,'D_LVIS', D_LVIS);
setappdata(gca,'D_06', D_06);

hax(2:4)=cheek_by_jowl(3, 1, [0.15 0.1 0.8 0.5]);
for kP=1:3
    axes(hax(kP+1));
    xy0=unique(round_to(D_06(kP).x(:), 25));
    [ddx, ddy]=meshgrid([-25:25:25]);
    ATM_mask=false(size(D_ATM.x));
    x0_ATM=round_to(D_ATM.x+1i*D_ATM.y, 25);
    for kk=1:length(ddx(:));
        ATM_mask=ATM_mask | ismember(x0_ATM, xy0+ddx(kk)+1i*ddy(kk));
    end
    LVIS_mask=false(size(D_LVIS.x));
    x0_LVIS=round_to(D_LVIS.x+1i*D_LVIS.y, 25);    
    for kk=1:length(ddx(:));
        LVIS_mask=LVIS_mask | ismember(x0_LVIS, xy0+ddx(kk)+1i*ddy(kk));
    end
    cla; hold on;
    
    if exist('D3','var')
        L=(real(D3(2*kP).x)-XL(1))/1000;
        plot(L, D3(2*kP).h_ph,'k.','markersize', 1)
    end
    
    L=(D_LVIS.x(LVIS_mask)-XL(1))/1000;
    plot_colored_points(L+1i*D_LVIS.z(LVIS_mask), D_LVIS.year(LVIS_mask), 2002:2020, [], this_cmap, true);
    L=(D_ATM.x(ATM_mask)-XL(1))/1000;
    plot_colored_points(L+1i*D_ATM.z(ATM_mask), D_ATM.year(ATM_mask), 2002:2020, [], this_cmap, true);
    L=(real(D_06(kP).x(:))-XL(1))/1000;
    hp_06=plot_colored_points( L+1i*D_06(kP).h_li(:), D_06(kP).year(:), 2002:2020, [], this_cmap, true);
    set(hp_06(ishandle(hp_06)),'markersize', 2)
    
end

set(hax(2:end),'yaxislocation','left');
ylabel(hax(3),'WGS-84 Elevation (m)');
set(hax(2:end),'tag', 'beam plots');

linkaxes(hax,'x');

h_bar=axes('position', [0.15 0.91 0.8 0.025]); hi=imagesc(2002:2020, [0 1], 2002:2020); 
colormap(h_bar, this_cmap);
set(h_bar,'ytick',[],'xaxislocation','top'); 
xlabel(h_bar,'year');

xlabel(hax(4),'x_{atc}, km');

set(hax(1),'visible','off')




 