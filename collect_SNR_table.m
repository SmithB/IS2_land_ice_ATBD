% load a first file to get the Hwin values
clear; load AT_fit_noise_test_June26_2018/test_file11.mat
% this gives us Hwin_vals and one BGR per file;

SNR_vals=-5:.1:5;

 
AT_fit_files=glob('AT_fit_noise_test_June26_2018/*.mat');
for k=1:length(AT_fit_files)
    temp=load(AT_fit_files{k},'BGR');
    BGR_vals(k)=temp.BGR;
end
[BGR_vals, ind]=sort(BGR_vals); 
AT_fit_files=AT_fit_files(ind);
SNR_F_table.P=zeros(length(SNR_vals), length(Hwin_vals), length(BGR_vals));


for kBG=1:length(AT_fit_files)
    load(AT_fit_files{kBG})
    for kH=1:length(Hwin_vals)
        if any(isfinite(out.P0(kH,:)))
            this_P0=median(out.P0(kH,isfinite(out.P0(kH,:))));
            for kSNR=1:length(SNR_vals)
                SNR_F_table.P(kSNR, kH, kBG)=this_P0*sum(isfinite(out.SNR(kH,:)) & out.SNR(kH,:)>SNR_vals(kSNR))/size(out.N,2);
            end
        end
    end
end

SNR_F_table.Hwin=Hwin_vals;
SNR_F_table.SNR=SNR_vals;
SNR_F_table.BGR=BGR_vals;
save SNR_table_june_2018 SNR_F_table

if false  % DON't do this until you are ready to overwrite the existing table
    out_file='SNR_F_table_June26_2018.h5';
    
    h5create(out_file,'/P_NoiseOnly', size(SNR_F_table.P),'datatype','double')
    h5write(out_file,'/P_NoiseOnly', SNR_F_table.P)
    
    h5create(out_file,'/SNR', size(SNR_F_table.P,1) ,'datatype','double')
    h5write(out_file,'/SNR', SNR_F_table.SNR)
    
    h5create(out_file,'/W_surface_window_initial', size(SNR_F_table.P,2),'datatype','double')
    h5write(out_file,'/W_surface_window_initial', SNR_F_table.Hwin)
    
    h5create(out_file,'/BGR',size(SNR_F_table.P,3),'datatype','double')
    h5write(out_file,'/BGR', SNR_F_table.BGR)
     
end



