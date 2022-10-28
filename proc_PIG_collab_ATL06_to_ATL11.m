 
IS_paths.ATL06='/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13C_NoFirn_NoDz/ATL06';
IS_paths.ATL11='/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13C_NoFirn_NoDz/ATL11/';

% load the ground tracks:
load data_files/PIG_groundtracks

% choose how many repeats to simulate:
Nreps=12;
% choose how many segments to include at each dh/dt point
Nsegs=6;
% how many times will we run this?  ATM, only one run exists
Nruns=1;

% choose whether to apply the delta_h function
APPLY_dh_function=false;
 
% define IS2 default parameters:
default_params.xt_scale=100;
default_params.beam_spacing=90;
default_params.time_zero=0;
default_params.sigma_geo=6.5;
default_params.t_scale=365;
t_launch=datenum('September 12 2018');

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
%for run_num=1:Nruns
run_num=1;
    run_dir=sprintf('%s/run_%d',  IS_paths.ATL11, run_num);
    if ~exist(run_dir,'dir')
        mkdir(run_dir);
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
 
            % now loop over repeats and load in the data  
            fprintf(1, 'ATL11_file is %s\n', ATL11_file);
            clear D3a; clear D3; clear D_ATL11a; count=0;
            for kR=1:Nreps
                % check if the ATL06 file exists
                ATL06_file=sprintf('%s/run_%d/rep_%d/Track_%d_D3.h5', IS_paths.ATL06, run_num, kR, TrackList(kT));
                if ~exist(ATL06_file,'file')
                    continue
                end
 
                % read the ATL06 file
                D3temp=read_ATL06_h5(ATL06_file, kP);
                if ~isfield(D3temp,'h_li'); continue; end
                if ~isfield(D3temp, 'atl06_quality_summary')
                    D3temp.ATL06_quality_summary=ATLAS_L3a_proc_ATBD('calc_ATL06_summary_flag',D3temp);
                end
                
                % remove the low-significance pairs
                good=any(D3temp.snr_significance < 0.05, 2);
                ff=fieldnames(D3temp);
                for kf=1:length(ff) 
                    D3temp.(ff{kf})=D3temp.(ff{kf})(good,:);
                end
                % set low-significance single-beam results to NaN;
                bad=D3temp.snr_significance > 0.05;
                h_fields={'h_li', 'h_li_sigma','h_mean','h_med'};
                for kf=1:length(h_fields)
                    D3temp.(h_fields{kf})(bad)=NaN;
                end
                
                % assign delta h
                if APPLY_dh_function
                    xx=ll2ps(D3temp.lat_ctr, D3temp.lon_ctr);
                    D3temp.delta_h=PIG_dh_function(real(xx), imag(xx), D3temp.time);
                    D3temp.h_li=D3temp.h_li+D3temp.delta_h;
                else
                    D3temp.delta_h=zeros(size(D3temp.delta_time));
                end
                
                % assign a missing parameter
                D3temp.ATL06_status=zeros(size(D3temp.x_atc));
                
                % interpolate the azimuth from the RGT
                D3temp.ref_azimuth=interp1(GroundTracks(TrackList(kT)).x_RGT, GroundTracks(TrackList(kT)).azimuth, D3temp.x_atc);
                
                % assign missing variable:
                pair=kP;
                D3temp.laser_beam=zeros(size(D3temp.h_li));
                D3temp.laser_beam(:,1)=2*pair-1;
                D3temp.laser_beam(:,2)=2*pair;
                % assign missing variable:
                D3temp.rep=ones(size(D3temp.h_li))*kR;
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
                 if ~ismember(f{kf},{'h_mean','h_med'})
                    D3.(f{kf})= cat(1, D3a.(f{kf}));
                 end            
            end
            [D3.x_ps_ctr, D3.y_ps_ctr]=ll2ps(D3.latitude, D3.longitude);
            D3.time=D3.delta_time;
            
            %  generate ATL11 structures
            clear D_ATL11a;
            % generate the fit centers:
            uX=min(D3.x_atc)+10*Nsegs:10*Nsegs:max(D3.x_atc)-10*Nsegs;
             
            D3.h_li_sigma=max(D3.sigma_h_mean, D3.fpb_med_corr_sigma);
            %assign cross-track slope errors
            D3.sigma_dh_fit_dy=repmat(sqrt(sum(D3.h_li_sigma.^2, 2))./diff(D3.y_atc, [], 2), [1, 2]);
            good=false(length(uX),1);
            % loop over along-track centers
            for kx=1:length(uX)
                %disp('WARNING: starting at 1000!!!!')
                % write ATL11-specific paramters into the _params_
                % structure
                this_param=params(1);
                this_param.x_atc_ctr=uX(kx);
                this_param.ref_pt_num=kx;
                this_param.max_degree_y=3;
                this_param.max_degree_x=4;
                
                % add uncorrected reflectance
                D3.n_noise=D3.w_surface_window_final.*D3.bckgrd/1.5e8;
                %D3.reflct_uncorr=(D3.n_fit_photons - D3.n_noise)./(D3.n_seg_pulses.*[3 12]);
                D3.reflct_uncorr=(D3.fpb_n_corr - D3.n_noise)./(D3.n_seg_pulses.*[3 12]);

                % select data around the current point
                D3b=index_struct(D3, abs(D3.x_atc-uX(kx)) < 20*Nsegs/2+2);
                
                % assign some extra fields
                fake_fields={'bsl_h','bslh_conf','cloud_flg','tide_ocean'};
                for kf=1:length(fake_fields)
                    D3b.(fake_fields{kf})=zeros(size(D3b.time));
                end
                D3b.sigma_g_h=zeros(size(D3b.time))+0.03;
                D3b.sigma_g_x=zeros(size(D3b.time))+6.5;
                D3b.sigma_g_y=zeros(size(D3b.time))+6.5;
                 
                % fit the data  HERE is the ATL11 processing 
                [D_ATL11a(kx), params]=ATL11_proc_ATBD_Nseg('calc_fit', D3b, Nreps, this_param);
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
            %top_groups={'corrected_h','reference_point','ref_surf','pass_stats','pass_quality','non_product'};
            top_groups={'corrected_h','reference_point','ref_surf','pass_stats','pass_quality'}
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
            write_ATL11_h5(ATL11_file, D_ATL11, D3, params);
          end
    end
%end
    