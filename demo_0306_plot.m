function demo_0306_plot(ATL03_file, x_atc_0, w)

%files=glob('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13B_NoFirn_NoDz_noDx/ATL03_subset/run_1/rep_*/Track_530-Pair_1_D2a.h5');
% note: for Track 530 Pair 1, rep 6 is great, rep 7 is so-so, 12 is pretty bad,
% smooth region;  2.8473e+07, +- 1 km
% demo_0306_plot(files{9}, 2.8473e+07, 2e3);
% 
% rough region:  2.8470e+07 +- 1.5 km




[D2a, PairData, params, TrackData]=read_ATLAS_h5_D2a(ATL03_file, true);
ATL06_file=strrep(ATL03_file,'ATL03','ATL06');
ATL06_file=strrep(ATL06_file, 'D2.h5','D3.h5');
[D3, dh_hist]=read_ATL06_h5(ATL06_file);
D3.atl06_quality_summary=ATLAS_L3a_proc_ATBD('calc_ATL06_summary_flag',D3);


figure; clf;
h=cheek_by_jowl(3,1, [0.15 0.07 0.7 0.86]);
colors={[.4 .4 1] [.9 0 0 ]};
for kB=1:2
    axes(h(kB));
    plot(D2a(kB).x_RGT(D2a(kB).SigNoise==1), D2a(kB).h(D2a(kB).SigNoise==1),'k.','markersize', 3);
    hold on;
    hs=plot_segs(D3.x_RGT(:, kB), D3.h_LI(:, kB), D3.dh_fit_dx(:, kB), 40, '-');
    set(hs,'color', colors{kB});
    
    %good=D3.SNR_significance(:, kB) < 0.02;
    good=D3.atl06_quality_summary(:, kB)==0;
    hs1=plot_segs(D3.x_RGT(good, kB), D3.h_LI(good, kB), D3.dh_fit_dx(good, kB), 40, '-');
    set(hs1(ishandle(hs1)),'linewidth', 3, 'color', colors{kB});  
    hs2=[plot_segs(D3.x_RGT(good, kB), D3.h_LI(good, kB)-D3.h_LI_sigma(good, kB), D3.dh_fit_dx(good, kB) , 40, '-');
        plot_segs(D3.x_RGT(good, kB), D3.h_LI(good, kB)+D3.h_LI_sigma(good, kB), D3.dh_fit_dx(good, kB), 40, '-')];
    set(hs2(ishandle(hs2)),'linewidth', 1, 'color', colors{kB});
end
 
D3.r_eff=(D3.n_fit_photons-D3.N_noise)/57; D3.r_eff(:,1)=D3.r_eff(:,1)/3; D3.r_eff(:,2)=D3.r_eff(:,2)/12;
axes(h(3)); hold on;
plot(D3.x_RGT(:,1), D3.r_eff(:,1),'r.');
plot(D3.x_RGT(:,2), D3.r_eff(:,2),'b.');
linkaxes(h,'x')

XY=ll2ps(D3.lat_ctr, D3.lon_ctr);

if exist('x_atc_0','var')
    set(h,'xlim', x_atc_0+[-w w]/2);
    XGL=get(gca,'xlim');
    els=D3.x_RGT(:,2) > XGL(:,1) & D3.x_RGT(:,2) < XGL(:,2);
    XR=range(real(XY(els,2)));
    YR=range(imag(XY(els,2)));
    W=max([diff(XR), diff(YR)])+500;
    DEM=read_geotif_xy('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_DEM/20140119_1453_10200100282DE200_102001002A716000-DEM_tr4x.tif', mean(XR)+[-W/2 W/2], YR+[-W/2 W/2]);
    DEM.z(DEM.z==0)=NaN;
    figure;
    hs=surf(DEM.x, DEM.y, DEM.z); shading interp; material dull; set(gca,'dataaspectratio', [1 1 .05]); colormap(jet*.6+.4); hl=light;
    
    hold on
    for kB=1:2
        e1=D2a(kB).x_RGT > XGL(:,1) & D2a(kB).x_RGT  < XGL(:,2);
        e1=e1 & D2a(kB).SigNoise==1;
        plot3(D2a(kB).xground(e1), D2a(kB).yground(e1), D2a(kB).zground(e1),'k.');
    end
end




