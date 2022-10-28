filenames=glob('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13B_NoFirn_NoDz/ATL06/run_1/rep_*/Track_462_D3.h5');

beams={'3l','3r'};
groups{1}='land_ice_height';
fields{1}={ 'atl06_quality_summary' 
         'delta_time' 
         'dh_fit_dx' 
         'dh_fit_dy' 
         'h_li' 
         'h_li_sigma' 
         'latitude' 
         'longitude' 
         'segment_id' 
         'sigma_geo_at' 
         'sigma_geo_h' 
         'sigma_geo_xt' };
groups{2}='ground_track';
fields{2}={'cycle','x_atc','y_atc'};
clear D6;
for kF=1:length(filenames)
    for kB=1:length(beams)
        for kg=1:length(groups)
            for kf=1:length(fields{kg})
                this_ds=sprintf('/gt%s/%s/%s/', beams{kB}, groups{kg}, fields{kg}{kf});
                D6(kF).(fields{kg}{kf})(:, kB)=h5read(filenames{kF}, this_ds);
            end
        end
    end
end

clear D6_sub;
x0=33046250.;
ff=fieldnames(D6);
for k=1:length(ff)
    D6_sub.(ff{k})=cat(1, D6.(ff{k}));
end
these= any(abs(D6_sub.x_atc(:,2)-x0)<=125, 2);
for k=1:length(ff)
    D6_sub.(ff{k})= D6_sub.(ff{k})(these,:);
end
 

