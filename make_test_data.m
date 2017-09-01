function varargout=make_test_data(varargin)

if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end

%----------------
function make_test_files

load WF_est
sigma0=1.6e-9*1.5e8;
Rough0=([.25 0.5, 1, 2, 4, 8])*sigma0;
Rsurf0=[0.0625 0.125 0.25 0.5 1];
BGR=1e7;
test_data_dir=['/home/ben/Dropbox/projects/IS2_ATBD/Combined_corr_test_data_Aug_30_2017'];
out_dir=test_data_dir;


% read in the SNR F table:
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields) 
    SNR_F_table.(fields{kf})=h5read('SNR_F_table.h5', ['/',fields{kf}]); 
end


GETLOG=true;

N_chan=[4 16];
N_per_white=N_chan*3/4;

if ~exist(test_data_dir,'dir')
    mkdir(test_data_dir);
    mkdir([test_data_dir,'/D2']);
    mkdir([test_data_dir,'/D3']);
end

[Rough,Rsurf]=meshgrid(Rough0, Rsurf0);
%Nsegs=10;
Nsegs=20;
N_pulses=Nsegs*58;
for k=1:numel(Rough)
    D2_file=sprintf('%s/D2/Rough=%3.2e_Rsurf=%3.2e_BGR=%3.2e.h5', test_data_dir, Rough(k), Rsurf(k), BGR);
    if exist(D2_file,'file')
        [D2, params]=read_test_file(D2_file);
        for beam=1:2; params(beam).WF=WF;end
        if GETLOG
            for beam=1:2
                params(beam).SAVELOG=true;
            end
        end
    else
        clear D2 params;
        for beam=1:2
            [D2(beam), params(beam)]=make_ATL03_data(N_pulses, N_chan(beam), Rough(k), N_per_white(beam)*Rsurf(k), WF, BGR, 40, 0.0125);
            %D2(beam)=index_struct(D2(beam), (D2(beam).h -D2(beam).z0)> -20 & (D2(beam).h-D2(beam).z0) < 30);  % random truncation...          
        end
         
        % assign left-out parameters
        for beam=1:2
            D2(beam).x_RGT=D2(beam).pulse_num*.7;
            D2(beam).y_RGT=zeros(size(D2(beam).x_RGT))+45*beam-90;
            D2(beam).time=D2(beam).pulse_num*1e-4;
            D2(beam).ph_class=zeros(size(D2(beam).h));
            XR=round_to(range(D2(beam).x_RGT), 10)+[-10 10];
            YR=round(range(D2(beam).h))+[-1 1];
            [count, xg, yg]=point_count_image(D2(beam).x_RGT, D2(beam).h, [XR(1) 5 XR(2)], [YR(1) 0.5 YR(2)]);  
            count=conv2(conv_corrected(full(count), gaussian([-6:6], 0, 2)), gaussian([-6:6]', 0, 1),'same');
            mask=conv2(double(count>1), ones(20,1),'same')>0;
            PCV=interp2(xg(:)', yg(:), count, D2(beam).x_RGT, D2(beam).h);
            maskval=interp2(xg(:)', yg(:), double(mask), D2(beam).x_RGT, D2(beam).h)==1;
            %             if (Rsurf(k) >=.25 && Rough(k) < 0.5 && beam ==1) | (Rsurf(k) >=.25 && Rough(k) <= 1 && beam ==2)
            %                 % if there's good signal strength, pick out the ground bin,
            %                 % OTW mark all as bad (use backup signal finder)
            %
            %                 D2(beam).ph_class(abs(D2(beam).h-0.5-D2(beam).z0) <10)=1;
            %                 D2(beam).ph_class(abs(D2(beam).h-0.5-D2(beam).z0)<2)=2;
            %             end
            if true
                D2(beam).ph_class(maskval==1)=1;
                D2(beam).ph_class(PCV>=1)=2;
            end
            %D2(beam).ph_class(abs(D2(beam).h)<2*sqrt(1^2+Rough(k)^2))=2;
            D2(beam).beam=zeros(size(D2(beam).h))+beam;
            D2(beam).track=zeros(size(D2(beam).h))+1;
            D2(beam).segment_number=floor(D2(beam).x_RGT/20)+1;
            % delete parameters that might cause confusion
            x_ps=D2(beam).x_RGT+ll2ps(-82, 116);
            [D2(beam).lat, D2(beam).lon]=ps2ll(x_ps);
        end
        for beam=1:2
            params(beam).N_channels=N_chan(beam);
            params(beam).roughness=Rough(k);
            params(beam).R_surf=Rsurf(k);
            if GETLOG
                params(beam).SAVELOG=true;
            end
        end
        
        for kB=1:2
            D2(kB)=index_struct(D2(kB), D2(kB).detected);
        end
        
        % delete problem parameters
        non_atl03_parameters={'x0','zground','xground','SigNoise','detected'};
        for k3=1:length(non_atl03_parameters)
            D2=rmfield(D2, non_atl03_parameters{k3});
        end    
        
        write_test_file(D2_file, D2, params);
    end
    %TEST_CODE% [D2, params]=read_test_file(D2_file);         for beam=1:2; params(beam).WF=WF;end
    if GETLOG
        [D3a, dh_hist, LOG]=ATLAS_L3a_proc_ATBD(D2, params, [], SNR_F_table);
    else
        [D3a, dh_hist]=ATLAS_L3a_proc_ATBD(D2, params, [], SNR_F_table);
    end
    clear D3;
    
    D3=D3a;
    
    D3=index_struct(D3, all(D3.first_seg_pulse > 80/.7 & D3.first_seg_pulse +58 < N_pulses - 80/.7, 2));
    D3_file=sprintf('%s/D3/Rough=%3.2e_Rsurf=%3.2e_BGR=%3.2e.h5', test_data_dir, Rough(k), Rsurf(k), BGR );
    write_struct_to_h5_file(D3_file, D3);
    if GETLOG
        LOGFILE=sprintf('%s/D3/Rough=%3.2e_Rsurf=%3.2e_BGR=%3.2e_LOG.mat', test_data_dir, Rough(k), Rsurf(k), BGR  );
        save(LOGFILE, 'LOG','dh_hist');
        LOG_SAVE=strrep(LOGFILE, '.mat','.h5');
        if exist(LOG_SAVE', 'file'); delete(LOG_SAVE); end
        dump_struct_to_h5(LOG, LOG_SAVE,'',{'seg%d_beam%d','it_%d'});

    end
end


%-------------------------------------------------------
function [D2, params]=make_ATL03_data(N_pulses, N_chan, roughness, N_per_pulse, WF, BGR, H_window, AT_slope)

if ~exist('H_window','var');
    H_window=4;
end

if exist('AT_slope','var');
    DEM=struct('x', -50:50:(50*ceil(N_pulses*.7/50+1)), 'y', [-100 100]');
    [xg,~]=meshgrid(DEM.x, DEM.y);
    DEM.z=xg*AT_slope;
else
    DEM=0;
end
x0=(1:N_pulses)*0.7;
ATM_xmit=ones(size(x0));
params=struct('roughness', roughness,   'sigma_x', 7.5, 'NoiseRate', BGR, 'H_window', H_window, 'WF', WF, 'c', 3e8, 'N_channels', N_chan, 'N_per_pulse', N_per_pulse,'t_dead', 3.2e-9);
D2=det_sim(DEM, x0, params, ATM_xmit);


%-------------------
function [f2, fp]=define_fields

f2={
    'segment_number'
    'pulse_num'
    'h'
    'z0'
    'BGR'
    'detected'
    'channel'
    'x_RGT'
    'y_RGT'
    'time'
    'ph_class'
    'beam'
    'track'
    'lat'
    'lon' };
fp={'H_window','N_channels','c','N_per_pulse','t_dead', 'sigma_x','roughness','R_surf'};

%-----------------------------------------------
function write_test_file(filename, D2,  params)

beam_name={'/weak/','/strong/'};
[f2, fp]=define_fields;
for k=1:length(f2)
    for beam=1:2
        if isfield(D2(beam),f2{k})
            write_h5_field(filename,  D2(beam).(f2{k}), ['/photon',beam_name{beam},f2{k}]);
        else
            warning('%s is not a field of D2', f2{k});
        end
    end
end
for k=1:length(fp)
    for beam=1:2
        if isfield(params(beam),fp{k})
            write_h5_field(filename,params(beam).(fp{k}), ['/params', beam_name{beam},fp{k}]);
        else
            warning('%s is not a field of params', fp{k});
        end
    end
end

%-------------------------------------------------------

function [D2, params]=read_test_file(filename)
beam_name={'/weak/','/strong/'};
[f2, fp]=define_fields;
for k=1:length(f2)
    for beam=1:2
        D2(beam).(f2{k})=h5read(filename, ['/photon', beam_name{beam},f2{k}]);
    end
end
for k=1:length(fp)
    for beam=1:2
        params(beam).(fp{k})=h5read(filename, ['/params', beam_name{beam},fp{k}]);
    end
end

%----------------------------------------
function write_struct_to_h5_file(filename, D3)
[thedir, thebase]=fileparts(filename);

if ~exist(thedir,'dir'); mkdir(thedir); end

if exist(filename,'file'); delete(filename); end
f_out=fieldnames(D3);
for k=1:length(f_out)
    write_h5_field(filename,D3.(f_out{k}), ['/',f_out{k}]);
end

%-------------------------------------------------------
function write_h5_field(filename, D, h5_field)
h5create(filename, h5_field, size(D),'datatype','double');
h5write(filename,h5_field, double(D));

