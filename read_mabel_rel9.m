function [D, GT]=read_mabel_rel9(filename, channel)

t0=(datenum(strrep(strrep(h5read(filename,'/ancillary_data/granule_start_utc'),'T',' '),'Z',' ')));


fields={'delta_time', 'ph_class','ph_h', 'ph_latitude','ph_longitude', 'ph_shot', 'ph_id'};

for kf=1:length(fields);
    DS=sprintf('/channel%03d/photon/%s', channel, fields{kf});
    D.(strrep(fields{kf},'ph_',''))=double(h5read(filename, DS));
end


D.time=(t0-datenum('jan 1 2000'))*24*3600+double(D.delta_time);

RE    = 6378206.4;
E2    = 0.006768658;

    
GT.latitude=h5read(filename,'/reference_track/rt_latitude');
GT.longitude=h5read(filename,'/reference_track/rt_longitude');

if max(D.latitude > 65); 
    D.x=gl_ll2ps(D.latitude, D.longitude);
    GT.x=gl_ll2ps(GT.latitude, GT.longitude);
    zone=1;
elseif min(D.latitude < -65);
    D.x=ll2ps(D.latitude, D.longitude); 
     GT.x=ll2ps(GT.latitude, GT.longitude);
    zone=-1;
else  % pick a UTM projection based on the first data point
    [N,E,zone]=ell2utm(D.latitude(1)*pi/180, D.longitude(1)*pi/180, RE, E2);
    lcm=deg2rad(zone*6-183);
    [N,E]=ell2utm(D.latitude*pi/180, D.longitude*pi/180, RE, E2, lcm);
    D.x=E+i*N;
    [N,E]=ell2utm(GT.latitude*pi/180, GT.longitude*pi/180, RE, E2, lcm);
    GT.x=E+i*N;
end

D.channel=channel*ones(size(D.latitude));
 
fields={'delta_time_end','delta_time_start','noise_rate'};
for kf=1:length(fields);
    DS=sprintf('/channel%03d/altimetry/%s', channel, fields{kf});
    Di.(strrep(fields{kf},'ph_',''))=double(h5read(filename, DS));
end
D.noise_rate=interp1((Di.delta_time_end+Di.delta_time_start)/2, Di.noise_rate, D.delta_time)*(7700/210);
D.noise_rate(~isfinite(D.noise_rate))=median(D.noise_rate(isfinite(D.noise_rate))) ;

 