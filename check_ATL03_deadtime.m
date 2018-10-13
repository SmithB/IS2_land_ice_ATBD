thefile='/Volumes/ice2/ben/scf/plateau_03/ATL03_2018_10_02_00100110_944_01.h5';

fields={'h_ph','ph_id_pulse','pce_mframe_cnt','ph_id_channel'};
beam='/gt1l';
clear D3;
for kf=1:length(fields);
    D3.(fields{kf})=double(h5read(thefile,[beam,'/heights/', fields{kf}]));
end


%[H,D3]=read_ATL03_photon_data('ATL03_2018_10_02_00100110_944_01.h5', 1, struct);
close all
% loop over channels in the beam:
uI=unique(D3(1).ph_id_channel);
for kCh=1:length(uI)
    % pull out the photons for the current channel
    D=index_struct(D3, D3.ph_id_channel==uI(kCh));
    % calculate the pulse number
    D.Pnum=double(D.ph_id_pulse)+200*(D.pce_mframe_cnt-min(D.pce_mframe_cnt));
    % sort by pulse number, then by (negative) height
    ph_scale=(max(D.h_ph)-min(D.h_ph))*1.05;
    [~, ind]=sort(D.Pnum+(max(D.h_ph)-D.h_ph)/ph_scale);
    D=index_struct(D, ind);
    
    % convert height to time
    dt=diff(D.h_ph/(-1.5e8));
    
    % select the delta-t values from the same pulse number
    good_dt=diff(D.Pnum)==0;
    dt=dt(good_dt);
    
    % make a histogram of the dt values.
    figure(double(uI(kCh)));
    histogram(dt, [-0.5e-8:1e-10:0.5e-8]);
end

