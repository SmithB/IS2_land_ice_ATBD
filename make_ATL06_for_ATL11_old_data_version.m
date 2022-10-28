in_top='/Volumes/ice1/ben/sdt/ATLxx_example/v13_hdf/PIG_ATL06_8.00MHz';
out_top='/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13/';
replace=true;
MAKE_NEW_ERROR_FILE=false;

firn_model_file=[out_top,'/zfirn_20km_10days_all_AA.mat'];
if ~exist('zfirn','var')
    load(firn_model_file);
end


firn_year_zero=1994.0;  % this for the firn model
launch_year=2018.75; 
delta_z_DEM=-8;

PairTrackCombos=load('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_pairtrack_list');
uTracks=unique(PairTrackCombos.track);
 
% setup the path defaults
IS_paths=IS_LI_paths('dh_clouds');
% Define the cloud fraction PDF
% optical thicknesses:
Pcloud.tau=0:4;
% probability of clouds at each thickness:
Pcloud.P=[24.6 37.5/2 37.5/2 37.9/2 37.9/2]/100;

% load the ground tracks:
load data_files/PIG_groundtracks

% choose how many repeats to simulate:
Nreps=12;
% choose how many times to run the simluation
N_runs=6;
for k_run=1:N_runs
    run_dir=sprintf('%s/run_%d', out_top, k_run); 
    if ~exist(run_dir,'dir'); mkdir(run_dir); end
    out_ATL06=[out_top,'/ATL06'];
    
    if ~exist(sprintf('%s/run_%d', out_ATL06, k_run),'dir'); mkdir(sprintf('%s/run_%d', out_ATL06, k_run)); end
    
    error_file{k_run}=sprintf('%s/ATL06/run_%d/applied_errors.mat', out_top, k_run);
    if ~exist(error_file{k_run},'file') || MAKE_NEW_ERROR_FILE     
        error_data.pointing.x=6.5*randn(Nreps, length(unique(PairTrackCombos.track)));
        error_data.pointing.y=6.5*randn(Nreps, length(unique(PairTrackCombos.track)));
        error_data.z=0.03*randn(Nreps, length(unique(PairTrackCombos.track)));
        error_data.track=unique(PairTrackCombos.track);
        error_data.rep=(1:Nreps)';
        save(error_file{k_run},'error_data');
    end
    tau_file{k_run}=sprintf('%s/ATL06/run_%d/tau_vals.mat', out_top, k_run);
    if ~exist(tau_file{k_run},'file')
        for kRep=1:Nreps
            for kTrack=1:length(uTracks)
                tau_ind=find(cumsum(Pcloud.P)>rand(1), 1, 'first');
                tau_vals.tau(kRep, kTrack)=Pcloud.tau(tau_ind);
            end  
        end
        tau_vals.track=uTracks(:)';
        tau_vals.Rep=(1:Nreps)';
        save(tau_file{k_run}, 'tau_vals');
    end
end
 
out_ATL06=[out_top,'/ATL06'];
for k_run=5:N_runs
    run_dir=sprintf('%s/run_%d', out_ATL06, k_run); if ~exist(run_dir,'dir'); mkdir(run_dir); end
    load(sprintf('%s/run_%d/tau_vals.mat', out_ATL06, k_run));
    load(sprintf('%s/applied_errors.mat', run_dir));
    for k_rep=1:Nreps
        out_rep_dir=sprintf('%s/rep_%d', run_dir, k_rep); if ~exist(out_rep_dir,'dir'); mkdir(out_rep_dir); end
         for k_TP=1:length(PairTrackCombos.track)
            track_ind=find(error_data.track==PairTrackCombos.track(k_TP));
            this_track=PairTrackCombos.track(k_TP);
            this_pair=PairTrackCombos.pair(k_TP);
            this_tau=tau_vals.tau(k_rep, track_ind);
            in_file=sprintf('%s/tau=%1.1f/rep_%d/Track_%d-Pair_%d_D3.h5', in_top, this_tau, k_rep, this_track, this_pair);
            out_file=sprintf('%s/Track_%d-Pair_%d_D3.h5', out_rep_dir, this_track, this_pair);
            if exist(out_file,'file') && replace==false; continue; end
            if ~exist(in_file,'file'); continue; end
            [D3, dh_hist]=read_ATL06_h5(in_file);
            x=ll2ps(D3.lat_ctr, D3.lon_ctr);
            % calculate zfirn 
            zfirn_interp=interpn(zfirn.y, zfirn.x, zfirn.year-min(zfirn.year)+firn_year_zero, zfirn.zs, imag(x), real(x),  D3.time/365.25+launch_year);
            % move the reported data locations
            x=x+error_data.pointing.x(k_rep,track_ind)+1i*error_data.pointing.y(k_rep,track_ind);
            
            % interpolate the azimuth from the RGT
            D3.azimuth=interp1(GroundTracks(this_track).x_RGT, GroundTracks(this_track).azimuth, D3.x_RGT);
   
            [D3.lat_ctr, D3.lon_ctr]=ps2ll(x);
            
            D3.delta_z_firn=zfirn_interp;
            D3.delta_z_error=error_data.z(k_rep, track_ind)*ones(size(D3.time));
            D3.delta_z_signal=PIG_dh_function(real(x), imag(x), D3.time);
            h_fields={'h_LI','h_mean','h_med'};
            
            for kf=1:length(h_fields)
                D3.(h_fields{kf})=D3.(h_fields{kf})+D3.delta_z_firn+D3.delta_z_error+D3.delta_z_signal+delta_z_DEM;
            end 
        
            D3.pointing_error_x=zeros(size(D3.h_LI))+error_data.pointing.x(k_rep,track_ind);
            D3.pointing_error_y=zeros(size(D3.h_LI))+error_data.pointing.y(k_rep,track_ind);
            write_ATL06_h5(D3, dh_hist, out_file);
        end
    end
end
    