function ATL06=ATL06_structure_gen(file_or_struct, dh_hist, ATL06_out_file, pair)
%%%%%%%%%%% MODIFIED From Kaitlin / Ben generated script
% (C) Nick Holschuh - U. of Washington - 2018 (Nick.Holschuh@gmail.com)
% This is designed to take synthetic matlab products and generated the
% expected hdf5 structure for processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The inputs are as follows:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% file_or_struct - Either a structure of size 3, containing the information
%                  for each beam pair, or the name of an h5 file to read in
% ATL06_out_file - string, containing the name of the output file with
%                  extension
% pair - vector containing any combination of [1 2 3], indicating the beam
%                  pair desired
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The outputs are as follows:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ATL06 - The final structure that is written to the h5 file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


beam_strings = [string('gt1l'),string('gt1r'),string('gt2l'),string('gt2r'),string('gt3l'),string('gt3r')];

if ~exist('pair','var')
    pair=[1 2 3]
end
beams=[2*pair-1; 2*pair]; beams=beams(:);


if isstruct(file_or_struct)
    ATL06_pull=file_or_struct;
else
    ATL06_file=file_or_struct;
    if strcmp(ATL06_file(end-2:end),'.h5')
        ATL06_pull = read_ATL06_h5(ATL06_in_file); % pull out the stuff in said simulated ATL06 file
    elseif strcmp(ATL06_file(end-3:end),'.mat')
        % this is a mat file, probably generated from an ATL03 by Ben.
        load(ATL06_file,'dh_hist','D_06');
        pair_opts = 1:length(D_06);
        beams = [];
        counter = 1;
        ATL06_pull = struct();
        for j = 1:length(pair_opts)
            if ~isempty(D_06{j})
                beams(counter) = j;
                if counter == 1
                    ATL06_pull= D_06{j};
                else
                    ATL06_pull(counter) = D_06{j};
                end
                counter = counter+1;
            end
        end 
    end
end

beams=beams(ismember(beams, cat(1, ATL06_pull.beam)));

beam_strings=beam_strings(beams);

for b = 1:length(beam_strings) % make a structure for each beam or set this to just equal one if you only want to extract one)
    
    bi = 2-mod(b,2);
    b2 = ceil(b/2);
    

    %ATL06_pull(b2).signal_selection_source(find(ATL06_pull(b2).signal_selection_source == 2)) = 1;


    if ~isfield(ATL06_pull(b2), 'ATL06_quality_summary')
        ATL06_pull(b2).ATL06_quality_summary=uint8(~(ATL06_pull(b2).h_robust_spread<1 &...
            ATL06_pull(b2).h_LI_sigma < 1 & ...
            ATL06_pull(b2).signal_selection_source <=1 & ...
            ATL06_pull(b2).SNR_significance < 0.02));

    end
    if size(ATL06_pull(b2).ATL06_quality_summary) ~= size(ATL06_pull(b2).h_LI_sigma)
        ATL06_pull(b2).ATL06_quality_summary=uint8(~(ATL06_pull(b2).h_robust_spread<1 &...
            ATL06_pull(b2).h_LI_sigma < 1 & ...
            ATL06_pull(b2).signal_selection_source <=1 & ...
            ATL06_pull(b2).SNR_significance < 0.02));        
    end
    
    if isfield(ATL06_pull(b2),'seg_count') == 0
        data.segment_quality.segment_id = round(ATL06_pull(b2).x_RGT(:,bi)/20);
    else
        data.segment_quality.segment_id = ATL06_pull(b2).seg_count(:,bi);
    end
    data.segment_quality.laser_beam = ATL06_pull(b2).beam(:,bi);
    data.segment_quality.reference_pt_lat =[];
    data.segment_quality.reference_pt_lon =[];
    data.segment_quality.delta_time = ATL06_pull(b2).time(:,bi);
    data.segment_quality.record_number = [];
    
    % N.B.  I've moved the signal_selection_status values up a level
    data.segment_quality.signal_selection_source = ATL06_pull(b2).signal_selection_source(:,bi);
    data.segment_quality.signal_selection_status_all = ATL06_pull(b2).signal_selection_status_all(:,bi);
    data.segment_quality.signal_selection_status_backup = ATL06_pull(b2).signal_selection_status_backup(:,bi);
    data.segment_quality.signal_selection_status_confident = ATL06_pull(b2).signal_selection_status_confident(:,bi);
    
  
    data.land_ice_segments.atl06_quality_summary = uint8(ATL06_pull(b2).ATL06_quality_summary(:,bi));
    data.land_ice_segments.delta_time=  ATL06_pull(b2).time(:,bi);
  
    data.land_ice_segments.h_li = ATL06_pull(b2).h_LI(:,bi);
    data.land_ice_segments.h_li_sigma = ATL06_pull(b2).h_LI_sigma(:,bi);
    data.land_ice_segments.latitude = ATL06_pull(b2).lat_ctr(:,bi);
    data.land_ice_segments.longitude = ATL06_pull(b2).lon_ctr(:,bi);
    if isfield(ATL06_pull(b2),'seg_count') == 0
        data.land_ice_segments.segment_id = round(ATL06_pull(b2).x_RGT(:,bi)/20);
    else
        data.land_ice_segments.segment_id = ATL06_pull(b2).seg_count(:,bi);
    end
    data.land_ice_segments.sigma_geo_h = sqrt((ATL06_pull(b2).sigma_geo_AT(:,bi).*ATL06_pull(b2).dh_fit_dx).^2 + ...
         (ATL06_pull(b2).sigma_geo_XT(:,bi).*ATL06_pull(b2).dh_fit_dy).^2 + ...
         0.03^2);
    
    data.bias_correction.fpb_mean_corr = ATL06_pull(b2).fpb_mean_corr(:,bi);
    data.bias_correction.fpb_mean_corr_sigma = ATL06_pull(b2).fpb_mean_corr_sigma(:,bi);
    data.bias_correction.fpb_med_corr_sigma = ATL06_pull(b2).fpb_med_corr_sigma(:,bi);
    data.bias_correction.fpb_med_corr = ATL06_pull(b2).fpb_med_corr(:,bi);
    data.bias_correction.fpb_n_corr = ATL06_pull(b2).fpb_N_corr(:,bi);
    data.bias_correction.med_r_fit = ATL06_pull(b2).med_r_fit(:,bi);
    data.bias_correction.tx_mean_corr = ATL06_pull(b2).TX_mean_corr(:,bi);
    data.bias_correction.tx_med_corr = ATL06_pull(b2).TX_med_corr(:,bi);
    data.bias_correction.tx_med_corr = ATL06_pull(b2).TX_med_corr(:,bi);
    
    data.dem.dem_flag = [];
    data.dem.dem_h = [];
    data.dem.geoid_h = [];
    
    data.fit_statistics.dh_fit_dx_sigma = ATL06_pull(b2).sigma_dh_fit_dx(:,bi);
    data.fit_statistics.h_expected_rms = ATL06_pull(b2).h_expected_rms(:,bi);
    data.fit_statistics.h_rms_misft = ATL06_pull(b2).h_rms(:,bi);
    data.fit_statistics.h_robust_sprd = ATL06_pull(b2).h_robust_spread(:,bi);
    data.fit_statistics.n_fit_photons = ATL06_pull(b2).n_fit_photons(:,bi);
    data.fit_statistics.n_seg_pulses = ATL06_pull(b2).N_seg_pulses(:,bi);
    data.fit_statistics.sigma_h_mean = ATL06_pull(b2).sigma_h_mean(:,bi);
    data.fit_statistics.signal_selection_source = ATL06_pull(b2).signal_selection_source(:,bi);
    data.fit_statistics.signal_selection_source_status = ATL06_pull(b2).signal_selection_status_all(:,bi);
    data.fit_statistics.snr = ATL06_pull(b2).SNR(:,bi);
    data.fit_statistics.snr_significance = ATL06_pull(b2).SNR_significance(:,bi);
    data.fit_statistics.w_surface_window_final = ATL06_pull(b2).w_surface_window_final(:,bi);
    data.fit_statistics.dh_fit_dx = ATL06_pull(b2).dh_fit_dx(:,bi);
    data.fit_statistics.dh_fit_dy = ATL06_pull(b2).dh_fit_dy(:,bi);
    data.fit_statistics.h_mean=ATL06_pull(b2).h_mean(:, bi);
    
    
    data.ground_track.cycle = ATL06_pull(b2).cycle(:,bi);
    if isfield(ATL06_pull(b2),'track') == 0
        data.ground_track.gt = ones(size(ATL06_pull(b2).cycle(:,bi)));
    else
        data.ground_track.gt = ATL06_pull(b2).track(:,bi);
    end
    data.ground_track.laser_beam = ATL06_pull(b2).beam(:,bi);
    data.ground_track.ref_azimuth = [];
    data.ground_track.ref_coelv = [];

    
        
    data.ground_track.rgt = [full(ATL06_pull(b2).x_RGT(:,bi)) ATL06_pull(b2).y_RGT(:,bi)]; 
    
    nan_inds = find(isnan(ATL06_pull(b2).lat_ctr(:,bi)));
    
    
    ATL06_pull(b2).lat_ctr(nan_inds,bi) = -90;
    ATL06_pull(b2).lon_ctr(nan_inds,bi) = 0;
    data.ground_track.seg_azimuth =  platte_carre_WGS84([ATL06_pull(b2).lat_ctr(:,bi), ATL06_pull(b2).lon_ctr(:,bi)], [mean(ATL06_pull(b2).lat_ctr(:,bi)), mean(ATL06_pull(b2).lon_ctr(:,bi))], 'track_azimuth'); 
    data.ground_track.seg_azimuth(nan_inds) = NaN;
    
    data.ground_track.x_atc = ATL06_pull(b2).x_RGT(:,bi);
    data.ground_track.y_atc = ATL06_pull(b2).y_RGT(:,bi);
    data.ground_track.sigma_geo_at = ATL06_pull(b2).sigma_geo_AT(:,bi);
    data.ground_track.sigma_geo_xt = ATL06_pull(b2).sigma_geo_XT(:,bi);

    
    
    data.geophysical.bckgrd = ATL06_pull(b2).BGR(:,bi);
    data.geophysical.bsnow_conf = [];
    data.geophysical.bsnow_h = [];
    data.geophysical.cloud_flg_asr = [];
    data.geophysical.cloud_flg_atm = [];
    data.geophysical.e_bckgrd = [];
    data.geophysical.r_eff = [];
    data.geophysical.solar_azimuth = [];
    data.geophysical.solar_elevation = [];
    
    data.geophysical.dac = [];
    data.geophysical.neutat_delay_total = [];
    data.geophysical.tide_earth = [];
    data.geophysical.tide_load = [];
    data.geophysical.tide_ocean = [];
    data.geophysical.tide_pole = [];
    
    % residual histogram is a different size than the rest of the fields
    if ~isempty(dh_hist{bi, b2}) && isstruct(dh_hist{bi, b2})
        data.residual_histogram.bckgrd_per_bin = dh_hist{bi, b2}.bckgrd_per_bin;
        data.residual_histogram.count = dh_hist{bi,b2}.count;
        data.residual_histogram.delta_time = (-4.995:0.01:4.995)/(299792458/2);
        data.residual_histogram.dh = -4.995:0.01:4.995;
        data.residual_histogram.lat_mean = dh_hist{bi,b2}.lat_mean;
        data.residual_histogram.lon_mean = dh_hist{bi,b2}.lon_mean;
        data.residual_histogram.pulse_count = dh_hist{bi,b2}.N_pulses;
        data.residual_histogram.segment_id_list = dh_hist{bi,b2}.segment_id_list;
        data.residual_histogram.x_atc_mean = dh_hist{bi,b2}.x_RGT_mean;
    else
        data.residual_histogram=[]
    end
    ATL06.(beam_strings{b})=data;
    
    sz=size(data.land_ice_segments.h_li);
    
    groups={'land_ice_segments','segment_quality','residual_histogram'}; 
    for kg=1:length(groups)
        if ~isstruct(data.(groups{kg}))
            continue
        end
        fields=fieldnames(data.(groups{kg}));
        for kf=1:length(fields)
            DS=sprintf('/%s/%s/%s', beam_strings{b}, groups{kg}, fields{kf});
            temp=data.(groups{kg}).(fields{kf});
            if isempty(temp); temp=NaN(sz); end
            try
                % we'll try to create the group, and if this fails,
                % we'll just write to it anyway
                if min(size(temp))==1
                    h5create(ATL06_out_file, DS, length(temp),'ChunkSize', [min(size(temp,1), 1024)], 'Datatype','double','Deflate', 9);
                else
                    h5create(ATL06_out_file, DS, size(temp),'ChunkSize', [min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
                end
            end
            h5write(ATL06_out_file, DS, temp);
        end
    end
    
    subgroups={'bias_correction','fit_statistics','geophysical', 'ground_track'};
    for kg=1:length(subgroups)
        fields=fieldnames(data.(subgroups{kg}));
        for kf=1:length(fields)
            DS=sprintf('/%s/land_ice_segments/%s/%s', beam_strings{b}, subgroups{kg}, fields{kf});
            temp=data.(subgroups{kg}).(fields{kf});
            if isempty(temp); temp=NaN(sz); end
            try
                % we'll try to create the group, and if this fails,
                % we'll just write to it anyway
                if min(size(temp))==1
                    h5create(ATL06_out_file, DS, length(temp),'ChunkSize', [min(size(temp,1), 1024)], 'Datatype','double','Deflate', 9);
                else                    
                    h5create(ATL06_out_file, DS, size(temp),'ChunkSize', [min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
                end
            end
            h5write(ATL06_out_file, DS, temp);
        end
    end

    clear data
end


