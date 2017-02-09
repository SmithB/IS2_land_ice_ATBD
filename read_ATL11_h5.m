function [D_ATL11, D_ATL06]=read_ATL11_h5(filename)


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


