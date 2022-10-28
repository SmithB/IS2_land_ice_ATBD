if ~exist('MOG','var')
    MOG=read_geotif('/Volumes/ice1/ben/MOG/MOG_125.tif', 3, 3);
end
figure(1); 
image_struct(MOG); 
hold on;
for kf=1:length(files) 
    for beam={'gt1','gt2','gt3'}
        D=read_ATL06_h5(files{kf}, beam{1}, 1); 
        hp=plot_colored_points(D.x+1i*D.y, D.h_li_sigma, [0.01:.01:1]);
        set(hp(ishandle(hp)),'tag', sprintf('%d %s', kf, beam{1})); 
    end
end
