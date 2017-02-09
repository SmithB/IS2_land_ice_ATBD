function [D3a, D, dh_hist]=proc_matlas_for_lori(filename)

if ~exist(filename,'file')
    error('%s does not exist', filename);
end

D=read_matlas(filename);
D=D.center(1);
D.pulse_num=double(D.shot_num);
 
Noise=h5read(filename, '/ancillary_data/matlas/al2a_atlas_noise_rate')*1000;

temp_h5='/var/tmp/matlas.h5';

if exist(temp_h5,'file')
    delete(temp_h5);
end
fields=fieldnames(D);
for kf=1:length(fields)
    h5create(temp_h5,['/ch00/',fields{kf}],size(D.time));
    h5write(temp_h5, ['/ch00/',fields{kf}], D.(fields{kf}));
end

 
% run the signal finder
fid=fopen('/var/tmp/matlas.lst','w');
fprintf(fid,'/var/tmp\nmatlas.h5\n011\nch00\n10none\n');
fclose(fid);
! idl -e 'mabel_signal_finding_driver, h5format=1, plot=0, debug=0, h5out=1, listfile="/var/tmp/matlas_new.lst"'

 
D.BGR=zeros(size(D.time))+Noise;

clear params
[params(1:2)]=struct('N_channels', 8,'sigma_pulse', 0.7e-9,'BGR', Noise,'c', 1.5e8,'sigma_x', 5, 't_dead', 3e-9);

D.ph_class=h5read('/var/tmp/matlas.h5V052015.0_sigparms.h5','//channelch00/photon/ph_class');
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


[D3a, dh_hist]=ATLAS_L3a_proc_ATBD(D, params);
