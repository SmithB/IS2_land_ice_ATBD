function [D_ATL11, D_ATL06, params]=read_ATL11_h5(filename)


II=h5info(filename);
for kG=1:length(II.Groups);
    theGroup=II.Groups(kG);
    GroupName=theGroup.Name;
    groupFieldName=strrep(GroupName,'/','');
    clear temp;
    for kD=1:length(theGroup.Datasets);
        theDataset=theGroup.Datasets(kD);
        dataName=theDataset.Name;
        dataFieldName=strrep(dataName,'/','');
        temp.(dataFieldName)=h5read(filename,[GroupName,'/', dataName]);
    end
    D_ATL11.(groupFieldName)=temp;    
end

D_ATL06=D_ATL11.ATL06;
D_ATL11=rmfield(D_ATL11, 'ATL06');

% read attributes:
fields={'poly_exp_x','poly_exp_y', 'time_zero'};
for kf=1:length(fields);
    params.(fields{kf})=h5readatt(filename,'/ref_surf/',fields{kf}); 
end 
