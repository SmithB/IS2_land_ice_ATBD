function [D3, Dhist]=read_ATL06_h5(h5_file, pair, hemisphere)

if ~exist('pair','var')  % we're reading an old-format file with one pair per file
    
    I=h5info(h5_file,'/');
    for kD=1:length(I.Datasets)
        D3.(I.Datasets(kD).Name)=h5read(h5_file,['/', I.Datasets(kD).Name]);
    end
    
    D3.seg_count=repmat((1:size(D3.x_RGT,1))', [1, 2]);
    
    LR_names={'L','R'};
    
    % this will work only if there's a dh_hist group
    try
        I=h5info(h5_file,'/hist_200m/L');
        for kD=1:length(I.Datasets);
            for kLR=1:2;
                Dhist(kLR).(I.Datasets(kD).Name)=h5read(h5_file,['/hist_200m/', LR_names{kLR},'/', I.Datasets(kD).Name]);
            end
        end
    catch
        Dhist=[];
    end
    
elseif ~ischar(pair)
    D3=struct;
    pair_name=sprintf('PT_%d', pair);
    top_groups={'ground_track', 'land_ice_height','fit_statistics','bias_correction'};
    
    for kg=1:length(top_groups)
        this_group=sprintf('/%s/%s', pair_name, top_groups{kg});
        try
            I=h5info(h5_file,this_group);
        catch
            continue
        end
        for kDS=1:length(I.Datasets)
            this_DS=I.Datasets(kDS).Name;
            D3.(I.Datasets(kDS).Name)=h5read(h5_file, [this_group,'/', this_DS]);
        end
        D3.bckgrd=h5read(h5_file,['/',pair_name,'/geophysical/bckgrd']);
    end
else
    if pair(end)=='l' || pair(end)=='r'
        beam={pair(end)};
        pair=pair(1:end-1);      
    else
        beam={'l','r'};
    end
    for kB=1:length(beam)
        %fprintf(1, '%s%s\n', pair, beam{kB});
        sub_group{1}='land_ice_height';
        datasets{1}={'h_li','h_li_sigma','latitude','longitude', 'segment_id','dh_fit_dx','dh_fit_dy','atl06_quality_summary','delta_time'};
        sub_group{2}='fit_statistics';
        datasets{2}={'dh_fit_dx_sigma','h_rms_misft','n_fit_photons','w_surface_window_final','snr_significance'};
        sub_group{3}='ground_track';
        datasets{3}={'x_atc','y_atc'};
        for kg=1:length(sub_group)
            for kd=1:length(datasets{kg})
                this_ds=datasets{kg}{kd};
                D3.(this_ds)(:, kB)=h5read(h5_file, sprintf('/%s%s/%s/%s', pair, beam{kB}, sub_group{kg}, this_ds));
            end
        end
        if exist('hemisphere', 'var')
            if hemisphere==1
                [D3.x, D3.y]=gl_ll2ps(D3.latitude, D3.longitude);
            else
                [D3.x, D3.y]=ll2ps(D3.latitude, D3.longitude);
            end
        end
    end          
end
