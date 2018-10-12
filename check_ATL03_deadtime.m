[H,D3]=read_ATL03_photon_data('ATL03_2018_10_02_00100110_944_01.h5', 1, struct);
close all
% loop over beams
for kB=1:2
    % loop over channels in the beam:
    uI=unique(D3(1).ph_id_channel);
    for kCh=1:length(uI)
        % pull out the photons for the current channel
        D=index_struct(D3(kB), D3(kB).ph_id_channel==uI(kCh));
        % calculate the pulse number
        D.Pnum=double(D.ph_id_pulse)+200*double(D.pce_mframe_cnt-min(D.pce_mframe_cnt));
        % sort by pulse number, then by (negative) height
        [~, ind]=sort(D.Pnum+(max(double(D.h_ph))+500-double(D.h_ph))/10000);
        D=index_struct(D, ind);
        
        % convert height to time
        dt=diff(double(D.h_ph)/(-1.5e8));
        
        % select the delta-t values from the same pulse number
        good_dt=diff(D.Pnum)==0;
        dt=dt(good_dt);
        
        % make a histogram of the dt values.
        figure(double(uI(kCh)));
        histogram(dt, [-0.5e-8:1e-10:0.5e-8]);
    end
end
