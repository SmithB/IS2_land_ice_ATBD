function D=read_h5_group(file, group, D)
if ~exist('D','var')
    D=struct();
end
II=h5info(file, ['/',group]);
for field={II.Datasets.Name}
    D.(field{1})=double(h5read(file, ['/',group,'/', field{1}]));
    try
        ndv=h5readatt(file,['/',group,'/', field{1}],'_FillValue');
        D.(field{1})(D.(field{1})==ndv)=NaN;
    catch
    end
end