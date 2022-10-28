
% input parameters:
version='v13_hdf';
apply_GF_error=false;
regenerate_y_offsets=false; 
N_reps=12;
base_BGR=8e6;
replace=true;


if false
    rep_list=1:N_reps;
    tau_vals=0:4;
    %tau_vals=[0 0.5 1 1.5 2 2.5 3 3.5 4];
    TP_list=1:32;
end


IS_paths=IS_LI_paths;
clipped_groundtrack_dir=[IS_paths.data_root,'/../Antarctica'];
global_groundtrack_dir=[IS_paths.data_root,'/Global'];
DEM_file=[IS_paths.data_root,'/PIG_DEM/20140119_1453_10200100282DE200_102001002A716000-DEM_tr4x.tif'];
DEM_HP_width=100; %width for filtering the DEM. The high-passed DEM is assumed to move at the surface velocity, the low-passed DEM is assumed static
T_track=91/1387*24*3600; %Orbital period, in seconds. Incrementing the track number by one increases the time by this number
GT_speed=0.7/1e-4; % use GT speed of 0.7 m /0.0001 s


if ~exist('DEM','var')
    DEM=getappdata(0,'DEM');
    DEM_LP_sub=getappdata(0,'DEM_LP_sub');
    if isempty(DEM_LP_sub);
        % read and filter the DEM
        DEM=read_geotif(DEM_file);
        DEM.z(DEM.z==0)=NaN;
        dx=DEM.x(2)-DEM.x(1);
        %     DEM_LP=read_geotif(strrep(DEM_file,'.tif','_LP100m.tif'));
        %     DEM.Z_lp=DEM_LP.z;
        %     clear DEM_LP
        %DEM.Z_lp=single(conv_corrected(double(DEM.z), gaussian(-128*dx:dx:128*dx, 0, DEM_HP_width)));
        %DEM.Z_hp=double(DEM.z)-DEM.Z_lp;
        % DEM.gx_lp=gradient(DEM.Z_lp, DEM.x, DEM.y);
        temp=read_geotif(strrep(DEM_file,'.tif','_LP100m_16x.tif'));
        DEM.xsub=temp.x(1:4:end);
        DEM.ysub=temp.y(1:4:end);
        DEM.zsub=temp.z(1:4:end, 1:4:end);
        DEM_LP_sub=struct('xsub', DEM.xsub,'ysub', DEM.ysub, 'zsub', DEM.zsub);
        Ng=500/diff(DEM.xsub(1:2));
        DEM_LP_sub.zsub=conv_corrected(DEM_LP_sub.zsub, gaussian(-ceil(3*Ng):ceil(3*Ng), 0, Ng), true);
        setappdata(0,'DEM', DEM);
        setappdata(0,'DEM_LP_sub', DEM_LP_sub);
    end
end

if false
    eval('make_groundtracks_for_ATL11_example');
else
    load([IS_paths.data_root,'/PIG_groundtracks'], 'GroundTracks');
end

load([IS_paths.data_root,'/PIG_PairTracks'])
load([IS_paths.data_root,'/PIG_LaserTracks'])

clear D;
f=fieldnames(LaserTracks);
for kf=1:length(f)
    D.(f{kf})=cat(1, LaserTracks.(f{kf}));
end

good=isfinite(interp2(DEM.xsub, DEM.ysub, DEM.zsub, real(D.xy), imag(D.xy)));
good=conv2(double(good), ones(5,1))~=0;
D=index_struct(D, good);

%load([data_root,'/GAUSSIAN_WF_est'])
load([IS_paths.data_root,'/WF_est'])

params_L=struct('N_per_pulse', 12, 't_dead', 3.2e-9, 'sigma_x', 7.2,'sigma_pulse', 1.6e-9,'c', 3e8, 'N_channels', 16, 'NoiseRate', base_BGR, 'H_window', 600, 'WF', WF,'DEBUG', false);
params_R=params_L; params_R.N_per_pulse=3; params_R.N_channels=4;

D.pair=ceil(D.beam/2);
 
u_TP=unique([D.track, D.pair], 'rows');
TP_list=1:length(u_TP);

interp_fields={'xy','t','x_RGT'};
copy_fields={'track','beam','pair'};
clear D2a_out params;
if regenerate_y_offsets
    y_offsets=45*randn(max(u_TP(:,1)),20);
else
    load([IS_paths.data_root,'/y_offsets.mat'])
end


 
%tau_vals=[0 1 2];
%tau_vals=[3];
%u_TP=u_TP(12,:);

if apply_GF_error
    GF_error_mag=min(80, params_L.H_window/4);
    GF_error_file=sprintf('GF_error_%2.0f.mat', GF_error_mag);
    if exist(GF_error_file,'file');
        load(GF_error_file);
    else
        GF_error_scale=ceil(5000/diff(DEM.xsub(1:2)));
        GF_error_dec=floor(GF_error_scale/8);
        % generate the GF_error for the different repeats and track/pair combos
        for k_rep=1:20
            for k_TP=1:length(u_TP)
                K_error=gaussian(-3*GF_error_scale:3*GF_error_scale, 0, GF_error_scale);
                temp=conv2_separable(randn(size(DEM.zsub)), K_error/sum(K_error));
                temp=GF_error_mag*temp./std(temp(isfinite(temp)));
                temp(temp>2*GF_error_mag-10)=2*GF_error_mag-10;
                temp(temp < -2*GF_error_mag+10)=-2*GF_error_mag+10;
                GF_error(k_rep, k_TP)=struct('x', DEM.xsub(1:GF_error_dec:end),'y', DEM.ysub(1:GF_error_dec:end),'z',  temp(1:GF_error_dec:end, 1:GF_error_dec:end));
            end
        end
        save(GF_error_file, '-v7.3', 'GF_error');
    end
end



for rep_ind=1:length(rep_list)
    k_rep=rep_list(rep_ind);    
    for k_tau=1:length(tau_vals)
    
         %y_offsets=randn(max(u_TP(:,1)),1);
        %rep_dir=sprintf('%s/ATLxx_example/v1.04/PIG_ATL03/rep_%d', data_root,  k_rep);
        rep_dir=sprintf('%s/%s/PIG_ATL03_8.00MHz/tau=%3.1f/rep_%d', IS_paths.data_root, version,  tau_vals(k_tau), k_rep);
        
        if ~exist(rep_dir,'dir')
            mkdir(rep_dir);
        end
        
        for k=1:length(TP_list)
            this_TP=u_TP(TP_list(k),:);
            % make code parallelizable
            out_file=sprintf('%s/Track_%d-Pair_%d_D2.h5',rep_dir, this_TP(1), this_TP(2));
%             if  ~replace
%                 if exist(out_file,'file'); continue; end
%                 status=lockfile_tool('lock', out_file);
%                 if status~=0
%                     continue
%                 end
%             end
            fprintf(1,'working on %s\n', out_file);
            
            if apply_GF_error
                this_GF_error=interp2(GF_error(k_rep, k).x(:)',GF_error(k_rep, k).y(:), GF_error(k_rep, k).z, DEM_LP_sub.xsub(:)', DEM_LP_sub.ysub(:));
                temp=conv_corrected(this_GF_error, gaussian(-30:30, 0, 16), true);
                this_GF_error(~isfinite(this_GF_error))=temp(~isfinite(this_GF_error));
                DEM.zsub=DEM_LP_sub.zsub+this_GF_error ;
            else
                DEM.zsub=DEM_LP_sub.zsub;
            end
            
            this_track=this_TP(1);
            this_y_offset=y_offsets(this_track, k_rep);
            els=D.track==this_TP(1) & D.pair==this_TP(2);
            D_TP=index_struct(D, els);
            uB=unique(D_TP.beam);
            clear D2a_out TrackData PairData
            for kB=1:length(uB)
                D1=index_struct(D_TP, D_TP.beam==uB(kB));
                si=min(real(D1.x_RGT)):0.7:max(real(D1.x_RGT)); si=si(:);
                ATM_xmit=ones(length(si),1)*exp(-tau_vals(k_tau));
                %ATM_xmit=ones(length(si)+57*2,1)+.5+2*conv(randn(length(si)+57*2,1), ones(57,1)/57, 'same');
                %ATM_xmit=ATM_xmit(57:(58+length(si)));
                
                D2=D1;
                for kf=1:length(interp_fields)
                    D2.(interp_fields{kf})=interp1(real(D1.x_RGT), D1.(interp_fields{kf}), si);
                end
                for kf=1:length(copy_fields)
                    D2.(copy_fields{kf})=zeros(size(si))+D1.(copy_fields{kf})(1);
                end
                
                t_hat=diff(D1.xy(1:2)); t_hat=t_hat/abs(t_hat);
                D2.xy=D2.xy+1i*t_hat*this_y_offset;
                D2.x_RGT=D2.x_RGT+1i*this_y_offset;
                
                XR=range(real(D2.xy(:)))+[-100 100];
                YR=range(imag(D2.xy(:)))+[-100 100];
                
                DEM_sub=subset_image(DEM, XR, YR);
                
                if mod(uB(kB), 2)==0
                    D_all=det_sim(DEM_sub, D2.xy, params_L, ATM_xmit);
                    these_params=params_L;
                else
                    D_all=det_sim(DEM_sub, D2.xy, params_R, ATM_xmit);
                    these_params=params_R;
                end
                %D2a=index_struct(D_all, D_all.detected);
                D2a=D_all;
                D2a.track=zeros(size(D2a.x0))+this_TP(1);
                D2a.beam=zeros(size(D2a.x0))+uB(kB);
                D2a.time=D2.t(D2a.pulse_num)/24/3600+91*(k_rep-1);
                D2a.x_RGT=D2.x_RGT(D2a.pulse_num);
                [D2a.lat, D2a.lon]=ps2ll(D2a.x0);
                D2a.azimuth=interp1(GroundTracks(this_TP(1)).x_RGT, GroundTracks(this_TP(1)).azimuth, real(D2a.x_RGT),'pchip');
                % Calculate the RPT coords --- no longer included
                % D2a.x_RPT=RGT_coords(D2a.x0, PairTracks(u_TP(k,1),ceil(uB(kB)/2)).xy, real(PairTracks(u_TP(k,1), ceil(uB(kB)/2)).x_RGT), 5000, 200);
                D2a_out(kB)=D2a;
                params(kB)=these_params;
                TrackData(kB)=GroundTracks(this_TP(1));
                PairData(kB)=PairTracks(this_TP(1),ceil(uB(kB)/2));
            end
            D2a=D2a_out;
            [out_dir, ~]=fileparts(out_file);
            cplx_fields={'x0', 'xground', 'x_RGT'};
            for kF=1:length(cplx_fields)
                y_field=strrep(cplx_fields{kF},'x','y');
                for kB=1:2
                    temp=D2a(kB).(cplx_fields{kF});
                    D2a(kB).(cplx_fields{kF})=real(temp);
                    D2a(kB).(y_field)=imag(temp);
                end
            end
            write_D2_HDF(struct('file', out_file,'D2a', D2a, 'params',params, 'TrackData', TrackData,'PairData' , PairData), out_dir)
            %save(out_file, 'D2a', 'params', 'TrackData','PairData');
            if ~replace
                lockfile_tool('unlock', out_file);
            end
        end
    end
end


    
    
    
    
