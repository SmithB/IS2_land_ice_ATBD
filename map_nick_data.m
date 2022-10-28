function D=map_nick_data(files)
 
gt={'gt1','gt2','gt3'};
for k=1:length(files)
    for k_gt=1:3
        D(k, k_gt)=read_ATL06_h5(files{k}, gt{k_gt}, 1); 
    end
end

temp=cat(1, D.x); XR=range(temp(:)); temp=cat(1, D.y); YR=range(temp(:));
DEM=read_geotif_xy('/Volumes/insar7/ben/ArcticDEM/FilledDEM_100m.tif', XR, YR);
  
for k=1:length(files)
    for k_gt=1:3
        if ~isempty(D(k, k_gt).x)
            D(k, k_gt).DEM=interp2(DEM.x, DEM.y, DEM.z, D(k, k_gt).x, D(k, k_gt).y);
        end
    end
end

return

ff=figure; hold on;
  
N0=[12 3];

for k=1:length(files)
    if ~ishandle(ff); break; end
    for k_gt=1:3
        for kB=1:2
            L0=D(k,k_gt).x_atc(:, kB); 
            if false
                z0=D(k,k_gt).h_li(:, kB)-D(k,k_gt).DEM(:, kB);           
                good=isfinite(L0) & isfinite(z0);
             
                [out, bin_ctr, bin_percentiles, N]=block_percentile_filter(L0(good), double(z0(good)), [0.86], 2000);
                [Ls, ii]=unique(L0);
                x0=interp1(Ls, D(k, k_gt).x(ii, kB)+1i*D(k, k_gt).y(ii, kB), bin_ctr);
                
                hp=plot_colored_points(x0,  bin_percentiles(:,1), 0:.125:2);
            end
            
            if true                
                z0=D(k,k_gt).n_fit_photons(:, kB);         

                good=isfinite(L0) & isfinite(z0);
                
                [out, bin_ctr, bin_percentiles, N]=block_percentile_filter(L0(good), double(z0(good)), [0.5], 1000);
                [Ls, ii]=unique(L0);
                x0=interp1(Ls, D(k, k_gt).x(ii, kB)+1i*D(k, k_gt).y(ii, kB), bin_ctr);
                hp=plot_colored_points(x0,  bin_percentiles(:,1)/(57*N0(kB)), 0:.05:0.5);
                
                set(hp(ishandle(hp)),'tag', sprintf('%d %d', k, k_gt),'buttondownfcn','s=get(gcbo,''tag'');disp(s); figure; plot_nick_profile(s, D)');
            end
        end
    end
end

% 
% for kf=1:size(D,1) 
%     for kp=1:3; D1=D(kf, kp); 
%         plot(D1.x(1:3:end), D1.y(1:3:end),'k.','markersize', 2,'tag', [files{kf},' ', num2str(kp)]); 
%     end 
% end