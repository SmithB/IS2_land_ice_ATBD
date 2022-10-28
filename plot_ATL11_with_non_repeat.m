function [D11, D6, D200, params, handles]=plot_ATL11_with_non_repeat(ATL11_file, DEM, years) 

t_launch=datenum('September 12 2018');
if ~exist('years','var');
    years=2000:2022;
end
[D11, D6, params]=read_ATL11_h5(ATL11_file);
[thedir, thebase]=fileparts(ATL11_file);
D200_file=[thedir,'/D200/', thebase,'.h5'];
II=h5info(D200_file);
for kD=1:length(II.Datasets)
    thename=II.Datasets(kD).Name;
    D200.(thename)=h5read(D200_file,['/',thename]); 
end
 
hax(1)=axes('position', [0.15 0.5, 0.7, 0.3]);

hold on;
[~, ind]=ismember(D11.non_repeat_data.ref_pt, D11.reference_point.seg_count);
[hp1]=plot_colored_points_inorder(D11.reference_point.x_ATC(ind)/1000+1i*D11.non_repeat_data.h_shapecorr, date2year(D11.non_repeat_data.time), years, [], jet(26)*.8, true);

good=D11.corrected_h.pass_quality_summary==0 & D11.pass_quality.min_SNR_significance<0.02;
[r, ~]=find(good); 
hp2=plot_colored_points_inorder(D11.reference_point.x_ATC(r)/1000+1i*D11.corrected_h.pass_h_shapecorr(good), date2year(t_launch+D11.corrected_h.mean_pass_time(good)), years,[],  jet(26)*.8, true);
set(hp2(ishandle(hp2)),'markersize', 1);
ylabel('WGS84 elevation, m');
set(gca,'xlim', range(D11.reference_point.x_ATC(r))/1000);

x11=ll2ps(D11.pass_stats.mean_pass_lat, D11.pass_stats.mean_pass_lon);
x200=ll2ps(D200.lat, D200.lon);
XR=range(real(x11(:)));
YR=range(imag(x11(:))); 
xlabel('Along-track distance, km');
first=find(isfinite(x11), 1,'first');
last=find(isfinite(x11), 1,'last');
if real(x11(last)) < real(x11(first))
    set(gca,'xdir','reverse');
end
caxis(range(years)); colormap(0.7*jet);
hb=colorbar('east'); ylabel(hb,'year');
hax(2)=axes('position', [0.15 0.15 0.7 0.25]);

%DEM=read_geotif_xy('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_DEM/20140119_1453_10200100282DE200_102001002A716000-DEM_tr4x_LP100m.tif', XR, YR);
gx =gradient(DEM.z, DEM.x, DEM.y);
imagesc(DEM.x, DEM.y,DEM.gx); 
colormap(gca,'gray'); caxis([-1 1]*0.1); axis xy equal;
hold on;

[hm1]=plot_colored_points_inorder(x200, date2year(D200.time), years, [], jet(22)*0.8, true);
[hm2]=plot_colored_points_inorder(x11(good), date2year(t_launch+D11.corrected_h.mean_pass_time(good)), years, [], jet(22)*0.8, true);
 set(hax(2),'visible','off');
handles=struct('hp1', hp1,'hp2', hp2,'hm1', hm1,'hm2', hm2);

if false
    years=2000:2022;
    t_launch=datenum('September 12 2018');

    ATL11_files=(glob('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13/ATL11//run_1/*.h5'));
    k=21; 
    %x0= 3.31067e+07
     
    %x0=3.31335e+07 % good sampling of OIB data:
    
     %x0=3.307572295194053e+07
         x0= 3.312360482667039e+07

     %x0= 33082520 ;% IS1 crossover
    figure(k); clf; 
    set(gcf,'position', [60 60 730 1000]); 
    [D11, D6, D200,~, handles]=plot_ATL11_with_non_repeat(ATL11_files{k}, DEM, years);
     
    [iNR, i11]=ismember(D11.non_repeat_data.ref_pt, D11.reference_point.seg_count);
    D11.non_repeat_data.x_ATC=NaN(size(D11.non_repeat_data.lat));
    D11.non_repeat_data.x_ATC(iNR)=D11.reference_point.x_ATC(i11);
    
    els_11=abs(D11.reference_point.x_ATC-x0)<30;
    els_NR=abs(D11.non_repeat_data.x_ATC-x0)<30;
    xy_ctr=mean(ll2ps(D11.corrected_h.ref_pt_lat(els_11), D11.corrected_h.ref_pt_lon(els_11)));
    els_200=abs(D200.x+1i*D200.y-xy_ctr)<125;
    set(gca,'xlim', get(gca,'xlim')+[-10 0])
    
    
    figure(100+k); clf;
    hax=cheek_by_jowl(2, 1, [0.15 0.15 0.7 0.7]);
    axes(hax(1));
    plot_colored_points(date2year(D200.time(els_200))+1i*D200.z(els_200), date2year(D200.time(els_200)), years, [], jet(length(years))*0.7, true);
    hold on;
    D6.delta_h=PIG_dh_function(D6.x_PS_ctr, D6.y_PS_ctr, D6.time*365.25);
    D6.year=date2year(D6.time+t_launch);
    els06= abs(D6.x_RGT-x0)<500;
    good06=D6.ATL06_quality_summary==0;
    %h_all=plot_colored_points(D6.year(els06)+1i*(D6.h_LI(els06)+D6.delta_h(els06)), D6.year(els06), years, [], jet(length(years))*0.7, true);
    h_all=plot_colored_points(D6.year(els06)+1i*(D6.h_LI(els06)+0*D6.delta_h(els06)), D6.year(els06), years, [], jet(length(years))*0.7, true);
    
    set(h_all(ishandle(h_all)),'markersize', 1);
    els06=els06 & D6.ATL06_quality_summary==0;
    %h_good=plot_colored_points(D6.year(els06)+1i*(D6.h_LI(els06)+D6.delta_h(els06)), D6.year(els06), years, [], jet(length(years))*0.7, true);
    h_good=plot_colored_points(D6.year(els06)+1i*(D6.h_LI(els06)+0*D6.delta_h(els06)), D6.year(els06), years, [], jet(length(years))*0.7, true);
    grid on;
    caxis(range(years)); hb(1)=colorbar('east'); colormap(jet*.8);
    
    axes(hax(2)); cla; hold on;
    y11=date2year(D11.corrected_h.mean_pass_time(:)+t_launch);
    good=D11.corrected_h.pass_quality_summary==0 & D11.pass_quality.min_SNR_significance<0.01;
    errorbar(date2year(D11.non_repeat_data.time(els_NR)), D11.non_repeat_data.z(els_NR), D11.non_repeat_data.sigma(els_NR),'k.');
    errorbar(y11(good&els_11), D11.corrected_h.pass_h_shapecorr(good & els_11), D11.corrected_h.pass_h_shapecorr_sigma(good & els_11),'k.')

    plot_colored_points(date2year(D11.non_repeat_data.time(els_NR))+1i*D11.non_repeat_data.z(els_NR), ...
        date2year(D11.non_repeat_data.time(els_NR)), years, [], jet(length(years))*0.7, true);
    good=D11.corrected_h.pass_quality_summary==0 & D11.pass_quality.min_SNR_significance<0.01;
    y11=date2year(D11.corrected_h.mean_pass_time(:)+t_launch);
    plot_colored_points(y11(good&els_11)+1i*D11.corrected_h.pass_h_shapecorr(good & els_11), y11(good&els_11),  years, [], jet(length(years))*0.7, true);
    caxis(range(years)); hb(2)=colorbar('east'); colormap(jet*.8)
    YL=get(hax(2),'ylim'); set(hax,'ylim', YL+[-20 20]);
    set(hax(2),'yaxislocation','left')
    xlabel(hax(2),'year');
    set(hax,'xlim', [2002 2025]); grid on;
    set(hax(1),'xticklabel', []); 
     for kk=1:2;
        pos=get(hb(kk),'position'); 
        set(hb(kk),'position', [pos(1)+pos(3), pos(2), pos(3)/2, pos(4)]); 
        ylabel(hax(kk),'WGS84 elevation');
    end

    
    
end

