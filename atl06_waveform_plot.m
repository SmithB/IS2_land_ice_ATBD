function [D_06, dh_hist, dhR]=atl06_waveform_plot(thefile)
fprintf(1,'working on %s', thefile);

if ~exist('D_06','var')
    [D_06, dh_hist]=read_ASAS_ATL06(thefile);
end
 


latR=range(D_06(2).latitude(isfinite(D_06(2).latitude)));
if max(latR) > 0
    hemisphere=1;
else
    hemisphere=-1;
end

for kP=1:3
    if hemisphere==1
        D_06(kP).x=gl_ll2ps(D_06(kP).latitude, D_06(kP).longitude);
    else
        D_06(kP).x=ll2ps(D_06(kP).latitude, D_06(kP).longitude);
    end
    
    browse_axes=findobj('tag','browse_axes');
    if ishandle(browse_axes)
        XL=get(browse_axes,'xlim');
        YL=get(browse_axes,'ylim');
        xx=D_06(kP).x;
        els=real(xx) > XL(1) & real(xx) < XL(2) & imag(xx) > YL(1) & imag(xx) < YL(2);
        els=any(els,2);
        D_06(kP)=index_struct(D_06(kP), els);
    end
end

 
n_chan=[4 16];

dh_hist=resample_residual_histogram(dh_hist);

for kB=1:numel(dh_hist)
    if hemisphere==1
        dh_hist(kB).x=gl_ll2ps(dh_hist(kB).lat_mean, dh_hist(kB).lon_mean);
    else
        dh_hist(kB).x=ll2ps(dh_hist(kB).lat_mean, dh_hist(kB).lon_mean);
    end
    browse_axes=findobj('tag','browse_axes');
    if ishandle(browse_axes)
        XL=get(browse_axes,'xlim');
        YL=get(browse_axes,'ylim');
        xx=dh_hist(kB).x;
        els=real(xx) > XL(1) & real(xx) < XL(2) & imag(xx) > YL(1) & imag(xx) < YL(2);
        for field = {'segment_id','lat_mean','lon_mean','x_atc_mean','x','pulse_count', 'bckgrd_per_bin'}
            dh_hist(kB).(field{1})=dh_hist(kB).(field{1})(els);
        end
        for field = {'segment_id_list','count'}
            dh_hist(kB).(field{1})=dh_hist(kB).(field{1})(:,els);
        end
    end
 
    
end



colors={[0.7 .2 0],[.0 .2 .8]};
% loop over pairs, make one figure per pair
for kP=1:3
    figure(kP+10);clf;
    
    hax(1)=axes('position', [0.1 0.6 0.8 0.3]);
    hax(2)=axes('position', [0.1 0.35 0.8 0.2]);
    hax(3)=axes('position', [0.1 0.1 0.8 0.2]);
    
    %qs=D_06(kP).h_robust_sprd>1 | D_06(kP).h_li_sigma > 1 | D_06(kP).signal_selection_source >1 | D_06(kP).snr_significance > 0.02;
    qs=D_06(kP).h_robust_sprd>1 | D_06(kP).h_li_sigma > 1  | D_06(kP).snr_significance > 0.02;
    
    D_06(kP).atl06_quality_summary=qs;    
    good_AT=ATL06_AT_filter(D_06(kP), 2);
    good_cloud=false(size(D_06(kP).h_li));
    cloud_rat=D_06(kP).n_fit_photons./D_06(kP).w_surface_window_final;
    good_cloud(cloud_rat(:,1)>1,1)=1;
    good_cloud(cloud_rat(:,2)>4,1)=1;
    
    axes(hax(1));
    hold on;
    title({strrep(thefile,'_','\_'), sprintf('pair %d', kP)});
    for kB=2:-1:1
        % plot of height vs x_ATC
        axes(hax(1));
        plot(D_06(kP).x_atc, D_06(kP).h_li(:, kB),'.','color', 'k','markersize', 3);
        good=D_06(kP).atl06_quality_summary(:, kB)==0 & good_AT(:, kB) & good_cloud(:, kB);
        plot(D_06(kP).x_atc(good, kB), D_06(kP).h_li(good, kB),'o','color', colors{kB});
        
        % plot of photons/segment
        axes(hax(2));
        hold on;
        plot(D_06(kP).x_atc, D_06(kP).n_fit_photons(:, kB)./(D_06(kP).n_seg_pulses*n_chan(kB)),'.','markersize', 1,'color', colors{kB});
        set(gca,'ylim', [0 0.8]);
    end   
    ylabel(hax(1),'h, m');
    ylabel(hax(2),{'signal strength,' 'Ph/pulse/pixel'})
    
    axes(hax(3));
    temp=(dh_hist(kB, kP).count);
    temp(~isfinite(temp))=0;
    Ps=conv2(temp, ones(40,1));
    good_cols=isfinite(dh_hist(kB, kP).x_atc_mean);
    imagesc(dh_hist(kB, kP).x_atc_mean(good_cols), dh_hist(1,1).dh(:,1), log10(Ps(:, good_cols)));
    setappdata(gcf,'D', dh_hist(1:2, kP));
    ylabel('\delta h, m');
    axis xy
    colorbar('east');
    colormap([0 0 0; jet(255)*.8+.2]);

    xlabel(hax(3),'x_{atc}');
    setappdata(gcf,'D',dh_hist(:, kP))
    set(findobj(gcf,'type','axes'), 'buttondownfcn', 'WF_buttondown');
    set(hax,'yaxislocation','left');

    linkaxes(hax,'x');
    D_scrubbed(kP)=D_06(kP);
    good=D_06(kP).atl06_quality_summary==0 & good_AT;
    D_scrubbed(kP).h_li(~good)=NaN;
end    
 

for kP=1:3;
    D_06(kP).latitude(D_06(kP).latitude==0)=NaN;
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
    MOS='/Volumes/ice1/ben/MOA/moa_2009_1km.tif';
end
II=read_geotif_xy(MOS, XC+[-HW HW], YR+[-HW HW]);
figure(4); clf; image_struct(II); colormap(gray); hold on; plot(xx,'r.'); caxis([14000 17000]);

 