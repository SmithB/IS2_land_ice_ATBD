function [D3, Dhist]=read_ATL06_h5(h5_file, pair)

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
    
else
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
end
             

