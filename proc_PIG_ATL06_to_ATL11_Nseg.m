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

 % get the list of ATL06 tracks
[s, out]=unix(sprintf('ls %s/*/*/*.h5', IS_paths.ATL06));
out=strsplit(deblank(out));
TrackNum=[];
for k=1:length(out)
    temp=regexp(out{k},'Track_(\d+)','tokens');
    if ~isempty(temp)
        TrackNum(k)=str2num(temp{1}{1});
    end
end
TrackList=unique(TrackNum(TrackNum~=0));

% iterate: repeat the simulation with different random clouds and position
% errors
for run_num=1:N_runs
    run_dir=sprintf('%s/run_%d',  IS_paths.ATL11, run_num);
    if ~exist(run_dir,'dir')
        mkdir(run_dir);
    end
    % if the track offsets exist, load them
    if exist([run_dir,'/track_offsets.mat'],'file')
        load([run_dir,'/track_offsets.mat'])
    elseif  exist([run_dir,'/track_offsets.mat'],'file')
        [s, out]=unix(['cp ' single_seg_dir,'/track_offsets.mat ', run_dir]);
        load([run_dir,'/track_offsets.mat']);
    else
        track_offsets=randn(1400,12)*default_params.sigma_geo*sqrt(2);
        save([run_dir,'/track_offsets.mat'],'track_offsets');
    end
    
    % load the file of optical thicknesses if it exists
    if exist([run_dir,'/TauList.mat'],'file')   
        load([run_dir,'/TauList.mat'],'TrackList','TauList')
    else
        % make the optical thickness
        clear TauList;
        for kT=1:length(TrackList)
            for kR=1:N_runs+2;
                TauList{kT, kR}=Pcloud.tau{find(rand(1)<cumsum(Pcloud.P), 1,'first')};
            end
        end
        save([run_dir,'/TauList.mat'],'TauList','TrackList');
    end
    
    % assign the default parameters
    clear params;  
    [params(1), params(2)]=deal(default_params);
    
    % loop over RGTs
    for kT=1:length(TrackList)
        % loop over pairs
        for kP=1:3
            % check if the output is under construction, skip if so
            ATL11_file=sprintf('%s/Track_%d-Pair_%d_D3b.h5', run_dir, TrackList(kT), kP);
            if lockfile_tool('lock', ATL11_file); continue; end
            if exist(ATL11_file,'file')
                lockfile_tool('unlock', ATL11_file)
                continue;
            end
            
            % now loop over repeats and load in the data  
            fprintf(1, 'ATL11_file is %s\n', ATL11_file);
            clear D3a; clear D3; clear D_ATL11a; count=0;
            for kR=1:Nreps
                % check if the ATL06 file exists
                ATL06_file=sprintf('%s/tau=%s/rep_%d/Track_%d-Pair_%d_D3.h5', IS_paths.ATL06, TauList{kT, kR}, kR, TrackList(kT), kP);
                if ~exist(ATL06_file,'file')
                    continue
                end
 
                % read the ATL06 file
                D3temp=read_ATL06_h5(ATL06_file);
                
                % assign delta h
                if APPLY_dh_function
                    xx=ll2ps(D3temp.lat_ctr, D3temp.lon_ctr);
                    D3temp.delta_h=PIG_dh_function(real(xx), imag(xx), D3temp.time);
                    D3temp.h_LI=D3temp.h_LI+D3temp.delta_h;
                else
                    D3temp.delta_h=zeros(size(D3temp.time));
                end
                
                % assign a missing parameter
                D3temp.ATL06_status=zeros(size(D3temp.x_RGT));
                
                % introduce the cross-track geolocation error
                this_track=unique(D3temp.track(isfinite(D3temp.track)));            
                D3temp=shift_ATL06(D3temp, GroundTracks(this_track), track_offsets(this_track, kR));
                
                % assign missing variable:
                pair=kP;
                D3temp.beam=zeros(size(D3temp.beam));
                D3temp.beam(:,1)=2*pair-1;
                D3temp.beam(:,2)=2*pair;
                % assign missing variable:
                D3temp.rep=ones(size(D3temp.z0))*kR;
                D3temp.cycle=D3temp.rep;
                
                count=count+1;
                % accumulate the output
                D3a(count)=D3temp;
            end
            % put together the fields from D3a into a single collection of
            % points
            clear D3;
            if ~exist('D3a','var'); lockfile_tool('unlock', ATL11_file); continue; end
            f=fieldnames(D3a);
            for kf=1:length(f)
                D3.(f{kf})=cat(1, D3a.(f{kf}));
            end
            
            %  generate ATL11 structures
            clear D_ATL11a;
            % generate the fit centers:
            uX=min(D3.x_RGT)+10*Nsegs:10*Nsegs:max(D3.x_RGT)-10*Nsegs;
            
             
            D3.h_LI_sigma=max(D3.sigma_h_mean, D3.fpb_med_corr_sigma);
            %assign cross-track slope errors
            D3.sigma_dh_fit_dy=repmat(sqrt(sum(D3.h_LI_sigma.^2, 2))./diff(D3.y_RGT, [], 2), [1, 2]);
            good=false(length(uX),1);
            % loop over along-track centers
            for kx=1:length(uX)
                % write ATL11-specific paramters into the _params_
                % structure
                this_param=params(1);
                this_param.x_RGT_ctr=uX(kx);
                this_param.ref_pt_num=kx;
                this_param.max_degree_y=3;
                this_param.max_degree_x=4;
                % select data around the current point
                D3b=index_struct(D3, abs(D3.x_RGT-uX(kx)) < 20*Nsegs/2+2);
                % fit the data  HERE is the ATL11 processing 
                D_ATL11a(kx)=ATL11_proc_ATBD_Nseg('calc_fit', D3b, Nreps, this_param);
                % check if the results are of the right shape.pwd
                
                if ~isempty(D_ATL11a(kx).corrected_h.ref_pt_lat)
                    good(kx)=true;
                else
                    good(kx)=false;
                end
            end
            D_ATL11a=D_ATL11a(good);
            
            %stick together the ATL11a structures to make ATL11
            clear D_ATL11
            % catenate together the valid entries:       
            top_groups={'corrected_h','reference_point','ref_surf','segment_stats','segment_geolocation'};
            for kg=1:length(top_groups);
                temp=[D_ATL11a.(top_groups{kg})];
                clear temp1;
                if isstruct(temp);
                    f=fieldnames(temp);
                    for kf=1:length(f);
                        temp1.(f{kf})=cat(1, temp.(f{kf}));
                    end
                    D_ATL11.(top_groups{kg})=temp1;
                end
            end
            % assign a segment count
            D_ATL11.reference_point.seg_count=find(good);
            
            % write out the results:
            write_ATL11_h5(ATL11_file, D_ATL11, D3);
            % unlock the file 
            lockfile_tool('unlock', ATL11_file)
        end
    end
end
    