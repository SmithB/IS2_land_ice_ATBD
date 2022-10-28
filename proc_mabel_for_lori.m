function [D3a, D, dh_hist]=proc_mabel_for_lori(filename, channel)
 

% read in the SNR F table:
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields);
    SNR_F_table.(fields{kf})=h5read('SNR_F_table.h5', ['/',fields{kf}]);
end

if ~exist(filename,'file')
    error('%s does not exist', filename);
end

D=read_mabel_rel9(filename, channel);
D.pulse_num=double(D.shot);

F12={'noise_rate','BGR'; 
    'class','ph_class'; 
};
 
for k=1:size(F12,1)
    D.(F12{k,2})=D.(F12{k,1});
    D=rmfield(D,F12{k,1});
end

 
clear params
[params(1:2)]=struct('N_channels', 8,'sigma_pulse', 0.7e-9,'BGR', median(D.BGR),'c', 1.5e8,'sigma_x', 5, 't_dead', 3e-9);

D.x_RGT=(D.delta_time-D.delta_time(1))*200;

D(2)=D;
 
zero_fields={'track','beam'};
for kB=1:2
    for kf=1:length(zero_fields)
        D(kB).(zero_fields{kf})=zeros(size(D(kB).h));
    end
end

D(1).y_RGT=-45+zeros(size(D(1).h));
D(2).y_RGT=45+zeros(size(D(2).h));


for kB=1:2
    D(kB).lat=D(kB).latitude;
    D(kB).lon=D(kB).longitude;
end


[D3a, dh_hist]=ATLAS_L3a_proc_ATBD(D, params, [], SNR_F_table);
ff=fieldnames(D3a); 
for kf=1:length(ff); 
    D3a.(ff{kf})=D3a.(ff{kf})(:,1);
end
D=D(1);


if false
    
     

    mabel_dir='/Volumes/ice1/ben/sdt/MABEL/';
    channels=[1 9 10 12];
    output_dir='/home/ben/Dropbox/temp/MABEL_for_Holly';
    files=dir([mabel_dir,'/*.h5']);
    for kF=1:length(files)
        figure(kF); hold on;
        for kC=1:length(channels)
            out_file=sprintf('%s/%s_ch%d.mat', output_dir,strrep(files(kF).name,'.h5',''), channels(kC))
             [D3a, D]=proc_mabel_for_lori([mabel_dir,'/', files(kF).name], channels(kC));
            plot(D3a.x_RGT, D3a.h_med,'o');
            save(out_file,'D3a');
        end
    end
    
end



