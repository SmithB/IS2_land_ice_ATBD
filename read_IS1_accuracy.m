fid=fopen('IS1_laser_accuracy.csv', 'r');
IS1_accuracy.campaign=strsplit(fgetl(fid),',');

temp0=fgetl(fid); 
fclose(fid);
temp0=strsplit(deblank(temp0),',');
for k=1:length(temp0)
    temp=str2num(strrep(temp0{k},'+-',' '));
    IS1_accuracy.mean_xy(k)=temp(1);
    IS1_accuracy.sigma_xy(k)=temp(2);
    IS1_accuracy.pointing_error_tot(k)=sqrt(sum(temp.^2));
end

IS1_accuracy=index_struct(IS1_accuracy, [1, 5:length(IS1_accuracy.campaign)]);