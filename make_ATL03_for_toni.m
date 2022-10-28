in_top='/Volumes/ice1/ben/sdt/ATLxx_example/v13_hdf/PIG_ATL03_8.00MHz';
out_top='/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13B_NoFirn_NoDz_noDx/';
replace=true;
MAKE_NEW_ERROR_FILE=false;
APPLY_FIRN=false;
APPLY_DH_FUNCTION=false;
APPLY_DELTA_X=false;

XR(1,:)=[-1594000, -1578000]; YR(1,:)=[-252000, -241000];
XR(2,:)=[-1560000 -1542000]; YR(2,:)=[ -244000, -235000];


firn_model_file=[out_top,'/zfirn_20km_10days_all_AA.mat'];
if ~exist('zfirn','var') && APPLY_FIRN
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
load([IS_paths.data_root,'/PIG_PairTracks'])

% choose how many repeats to simulate:
Nreps=12;
% choose how many times to run the simluation
N_runs=1;
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
out_ATL03=[out_top,'/ATL03_subset'];
for k_run=1:N_runs
    run_dir=sprintf('%s/run_%d', out_ATL03, k_run); if ~exist(run_dir,'dir'); mkdir(run_dir); end
    load(sprintf('%s/run_%d/tau_vals.mat', out_ATL06, k_run));
    load(sprintf('%s/run_%d/applied_errors.mat', out_ATL06, k_run));
    for k_rep=1:Nreps
        out_rep_dir=sprintf('%s/rep_%d', run_dir, k_rep); if ~exist(out_rep_dir,'dir'); mkdir(out_rep_dir); end
        Tracks=unique(PairTrackCombos.track);
        
        for kT=1:length(Tracks)
            this_track=Tracks(kT);
            track_ind=find(tau_vals.track==this_track);
            out_file=sprintf('%s/Track_%d_D2a.h5', out_rep_dir, this_track);
            this_tau=tau_vals.tau(k_rep, track_ind);
             
            for this_pair=1:3
                in_file=sprintf('%s/tau=%1.1f/rep_%d/Track_%d-Pair_%d_D2.h5', in_top, this_tau, k_rep, this_track, this_pair);
                
                if ~exist(in_file,'file'); continue; end;
                out_file=sprintf('%s/Track_%d-Pair_%d_D2a.h5', out_rep_dir, this_track, this_pair);
                if exist(out_file,'file') && replace==false; continue; end
 
                %D2=read_ATLAS_h5_D2a(in_file);
                [D2, PairData, params, TrackData]=read_ATLAS_h5_D2a(in_file);
                for kB=1:2
                    x=ll2ps(D2(kB).lat, D2(kB).lon);
                    in_box=false(size(x));
                    for k_reg=1:2
                        in_box=in_box | (real(x) > XR(k_reg,1) & real(x) < XR(k_reg,2) & imag(x) > YR(k_reg,1) & imag(x) < YR(k_reg,2));
                    end
                    D2(kB)=index_struct(D2(kB), in_box);
                    x=x(in_box);
                    if isempty(D2(kB).x0); continue; end
                    
                    % calculate zfirn
                    if APPLY_FIRN
                        zfirn_interp=interpn(zfirn.y, zfirn.x, zfirn.year-min(zfirn.year)+firn_year_zero, zfirn.zs, imag(x), real(x),  D2(kB).time/365.25+launch_year);
                    else
                        zfirn_interp=zeros(size(x));
                    end
                    
                    if APPLY_DELTA_X
                        % move the reported data locations
                        x=x+error_data.pointing.x(k_rep,track_ind)+1i*error_data.pointing.y(k_rep,track_ind);
                    end
                    
                    % interpolate the azimuth from the RGT
                    D2(kB).azimuth=interp1(GroundTracks(this_track).x_RGT, GroundTracks(this_track).azimuth, D2(kB).x_RGT);
                    
                    [D2(kB).lat, D2(kB).lon]=ps2ll(x);
                    
                    D2(kB).delta_z_firn=zfirn_interp;
                    D2(kB).delta_z_error=error_data.z(k_rep, track_ind)*ones(size(D2(kB).time));
                    if APPLY_DH_FUNCTION
                        D2(kB).delta_z_signal=0*PIG_dh_function(real(x), imag(x), D2(kB).time);
                    else
                        D2(kB).delta_z_signal=zeros(size(x));
                    end
                    h_fields={'h'};
                    
                    for kf=1:length(h_fields)
                        D2(kB).(h_fields{kf})=D2(kB).(h_fields{kf})+D2(kB).delta_z_firn+D2(kB).delta_z_error+D2(kB).delta_z_signal+delta_z_DEM;
                    end
                    
                    D2(kB).pointing_error_x=zeros(size(D2(kB).h))+error_data.pointing.x(k_rep,track_ind);
                    D2(kB).pointing_error_y=zeros(size(D2(kB).h))+error_data.pointing.y(k_rep,track_ind);
                end
     
                if ~isempty(D2(1).x0) || ~isempty(D2(2).x0)
                    write_D2_HDF(struct('file', out_file,'D2a', D2, 'params',params, 'TrackData', TrackData,'PairData' , PairData), out_rep_dir)
                end
            end
        end
    end
end
    