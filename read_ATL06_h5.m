function [D3, Dhist]=read_ATL06_h5(h5_file)


I=h5info(h5_file,'/');
for kD=1:length(I.Datasets);
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


