function ATL06=ATL06_structure_gen_original(ATL06, file_or_struct, dh_hist,  pair)

beam_strings = [string('gt1l'),string('gt1r'),string('gt2l'),string('gt2r'),string('gt3l'),string('gt3r')];

if isstruct(file_or_struct)
    ATL06_pull=file_or_struct;
else
    ATL06_file=file_or_struct;
    [ATL06_pull, dh_hist] = read_ATL06_h5(ATL06_in_file); % pull out the stuff in said simulated ATL06 file
    %temp=regexp(ATL06_in_file,'Pair_(\d)_D','tokens');
    %pair=str2double(temp{1}{1});   
end
 
  
beam_strings=beam_strings(2*(pair-1)+[1:2]);
if ~isfield(ATL06_pull, 'ATL06_quality_summary')
    ATL06_pull.atl06_quality_summary=uint8(ATLAS_L3a_proc_ATBD('calc_ATL06_summary_flag',ATL06_pull));
end

good=all(isfinite(ATL06_pull.h_mean),2);

for b = 1:2 % make a structure for each beam or set this to just equal one if you only want to extract one)
     
    data.segment_quality.segment_id = ATL06_pull.segment_id(:, b);
    data.segment_quality.laser_beam = ATL06_pull.beam(:, b);
    data.segment_quality.reference_pt_lat =ATL06_pull.lat_ctr(:, b);
    data.segment_quality.reference_pt_lon =ATL06_pull.lon_ctr(:, b);
    data.segment_quality.delta_time = ATL06_pull.time(:, b);
    data.segment_quality.record_number = ATL06_pull.segment_id(:, b);
    data.segment_quality.signal_selection_source = ATL06_pull.signal_selection_source(:, b);
    data.segment_quality.signal_selection_status.signal_selection_status_all = ATL06_pull.signal_selection_status_all(:, b);
    data.segment_quality.signal_selection_status.signal_selection_status_backup = ATL06_pull.signal_selection_status_backup(:, b);
    data.segment_quality.signal_selection_status.signal_selection_status_confident = ATL06_pull.signal_selection_status_confident(:, b);
    
     
    data.land_ice_height.atl06_quality_summary = ATL06_pull.atl06_quality_summary(good, b);
    data.land_ice_height.delta_time=  ATL06_pull.time(good, b);
    data.land_ice_height.dh_fit_dx = ATL06_pull.dh_fit_dx(good, b);
    data.land_ice_height.dh_fit_dy = ATL06_pull.dh_fit_dy(good, b);
    data.land_ice_height.h_li = ATL06_pull.h_LI(good, b);
    data.land_ice_height.h_li_sigma = ATL06_pull.h_LI_sigma(good, b);
    data.land_ice_height.latitude = ATL06_pull.lat_ctr(good, b);
    data.land_ice_height.longitude = ATL06_pull.lon_ctr(good, b);
    data.land_ice_height.segment_id = ATL06_pull.segment_id(good, b);
    data.land_ice_height.sigma_geo_at = ATL06_pull.sigma_geo_AT(good, b);
    data.land_ice_height.sigma_geo_xt = ATL06_pull.sigma_geo_XT(good, b);
    data.land_ice_height.sigma_geo_h = 0.03+zeros(size(ATL06_pull.sigma_geo_AT(good, b)));
    
    blank=NaN(size(ATL06_pull.h_LI(good, b)));
    data.land_ice_height.bias_correction.fpb_mean_corr = ATL06_pull.fpb_mean_corr(good, b);
    data.land_ice_height.bias_correction.fpb_mean_corr_sigma = ATL06_pull.fpb_mean_corr_sigma(good, b);
    data.land_ice_height.bias_correction.fpb_med_corr_sigma = ATL06_pull.fpb_med_corr_sigma(good, b);
    data.land_ice_height.bias_correction.fpb_med_corr = ATL06_pull.fpb_med_corr(good, b);
    data.land_ice_height.bias_correction.fpb_n_corr = blank;
    data.land_ice_height.bias_correction.med_r_fit = ATL06_pull.med_r_fit(good, b);
    data.land_ice_height.bias_correction.tx_mean_corr = blank;
    data.land_ice_height.bias_correction.tx_med_corr = blank;
    data.land_ice_height.bias_correction.tx_mean_corr = ATL06_pull.TX_mean_corr(good, b);
    data.land_ice_height.bias_correction.tx_med_corr = ATL06_pull.TX_med_corr(good, b);
    data.land_ice_height.bias_correction.tx_med_corr = ATL06_pull.TX_med_corr(good, b);
    
    data.land_ice_height.dem.dem_flag = blank;
    data.land_ice_height.dem.dem_h = blank;
    data.land_ice_height.dem.geoid_h = blank;
    
    data.land_ice_height.fit_statistics.dh_fit_dx_sigma = blank;
    data.land_ice_height.fit_statistics.h_expected_rms = ATL06_pull.h_expected_rms(good, b);
    data.land_ice_height.fit_statistics.h_rms_misft = ATL06_pull.h_rms(good, b);
    data.land_ice_height.fit_statistics.h_robust_spread = ATL06_pull.h_robust_spread(good, b);
    data.land_ice_height.fit_statistics.n_fit_photons = ATL06_pull.n_fit_photons(good, b);
    data.land_ice_height.fit_statistics.n_seg_pulses = ATL06_pull.N_seg_pulses(good, b);
    data.land_ice_height.fit_statistics.sigma_h_mean = ATL06_pull.sigma_h_mean(good, b);
    data.land_ice_height.fit_statistics.signal_selection_source = ATL06_pull.signal_selection_source(good, b);
    data.land_ice_height.fit_statistics.signal_selection_source_status = ATL06_pull.signal_selection_status_all(good, b);
    data.land_ice_height.fit_statistics.snr = ATL06_pull.SNR(good, b);
    data.land_ice_height.fit_statistics.snr_significance = ATL06_pull.SNR_significance(good, b);
    data.land_ice_height.fit_statistics.w_surface_window_final = ATL06_pull.w_surface_window_final(good, b);
    
    data.land_ice_height.ground_track.cycle = ATL06_pull.cycle(good, b);
    data.land_ice_height.ground_track.gt = ATL06_pull.GT(good, b);
    data.land_ice_height.ground_track.laser_beam = ATL06_pull.beam(good, b);
    data.land_ice_height.ground_track.ref_azimuth = blank;
    data.land_ice_height.ground_track.ref_coelv = blank;
    data.land_ice_height.ground_track.rgt = blank;
    data.land_ice_height.ground_track.seg_azimuth = blank;
    data.land_ice_height.ground_track.x_atc = ATL06_pull.x_RGT(good, b);
    data.land_ice_height.ground_track.x_atc = ATL06_pull.x_RGT(good, b);
    data.land_ice_height.ground_track.y_atc = ATL06_pull.y_RGT(good, b);
    
    data.land_ice_height.sun_clouds.bckgrd = ATL06_pull.BGR(good, b);
    data.land_ice_height.sun_clouds.bsnow_conf = blank;
    data.land_ice_height.sun_clouds.bsnow_h = blank;
    data.land_ice_height.sun_clouds.cloud_flg_asr = blank;
    data.land_ice_height.sun_clouds.cloud_flg_atm = blank;
    data.land_ice_height.sun_clouds.e_bckgrd = blank;
    data.land_ice_height.sun_clouds.r_eff = blank;
    data.land_ice_height.sun_clouds.solar_azimuth = blank;
    data.land_ice_height.sun_clouds.solar_elevation = blank;
    
    data.land_ice_height.tides_atmosphere.dac = blank;
    data.land_ice_height.tides_atmosphere.neutat_total_delay = blank;
    data.land_ice_height.tides_atmosphere.tide_earth = blank;
    data.land_ice_height.tides_atmosphere.tide_load = blank;
    data.land_ice_height.tides_atmosphere.tide_ocean = blank;
    data.land_ice_height.tides_atmosphere.tide_pole = blank;
    
    blank=NaN(size(dh_hist(b).N_pulses));
    data.residual_histogram.bkgrd_expected = blank;
    data.residual_histogram.count =dh_hist(b).count;
    data.residual_histogram.delta_time = blank;
    data.residual_histogram.lat_mean =  dh_hist(b).lat_mean;
    data.residual_histogram.lon_mean = dh_hist(b).lon_mean;
    data.residual_histogram.pulse_count = dh_hist(b).N_pulses;
    data.residual_histogram.segment_id_list = NaN(length(dh_hist(b).N_pulses), 10);% reshape(dh_hist(b).segment_id_list, 10, length(dh_hist(b).segment_id_list)/10);
    data.residual_histogram.x_atc_mean =  dh_hist(b).x_RGT_mean;
    
    beam_assign = beam_strings{b};
    ATL06.(beam_assign) = data; 
end




