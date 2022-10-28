% setup the path defaults
IS_paths=IS_LI_paths('dh_clouds');

% Define the cloud fraction PDF
% optical thicknesses:
Pcloud.tau={'0.0' '1.0' '2.0' '3.0' '4.0'};
% probability of clouds at each thickness:
Pcloud.P=[24.6 37.5/2 37.5/2 37.9/2 37.9/2]/100;

% load the ground tracks:
load data_files/PIG_groundtracks

% choose how many repeats to simulate:
Nreps=12;
% choose how many times to run the simluation
N_runs=10;
% choose how many segments to include at each dh/dt point
Nsegs=6;

% choose whether to apply the delta_h function
APPLY_dh_function=true;

% name these conditions
run_type='dh_clouds';


% define IS2 default parameters:
default_params.xt_scale=100;
default_params.beam_spacing=90;
default_params.time_zero=0;
default_params.sigma_geo=6.5;
default_params.t_scale=365;

PairTrack_list=load('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_pairtrack_list');
 
TrackList=unique(PairTrack_list.pair);

% iterate: repeat the simulation with different random clouds and position
% errors
for run_num=1:N_runs
    run_dir=sprintf('%s/run_%d',  IS_paths.ATL11, run_num);
    if ~exist(run_dir,'dir')
        mkdir(run_dir);
    end
    % make the offset files
    
    track_offsets=randn(1400,12)*default_params.sigma_geo*sqrt(2);
    save([run_dir,'/track_offsets.mat'],'track_offsets');
    
    
    % load the file of optical thicknesses if it exists
    
    % make the optical thickness
    clear TauList;
    for kT=1:length(TrackList)
        for kR=1:Nreps+2
            TauList{kT, kR}=Pcloud.tau{find(rand(1)<cumsum(Pcloud.P), 1,'first')};
        end
    end
    save([run_dir,'/TauList.mat'],'TauList','TrackList');
end