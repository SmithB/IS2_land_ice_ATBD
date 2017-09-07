function varargout=make_fpb_test_data(varargin)

if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end

%----------------
function make_test_files

c2=1.5e8;
sigma0=0.68e-9*1.5e8;

WF.t=[-15e-9:1e-12:15e-9];
WF.p=gaussian(WF.t, 0, sigma0/c2);

Rough0=([ 0 1, 2, 4, 8])*sigma0;
Rsurf0=[0.25 0.5  0.75 1];
BGR=0;%1e7;
test_data_dir=['/home/ben/Dropbox/projects/IS2_ATBD/deadtime_correction_data/Sep_7_2017'];
out_dir=test_data_dir;
 
N_chan=[4 16];
N_per_white=N_chan*3/4;

if ~exist(test_data_dir,'dir')
    mkdir(test_data_dir);
end
if ~exist([test_data_dir,'/D2'],'dir')
    mkdir([test_data_dir,'/D2']);
end

[Rough,Rsurf]=meshgrid(Rough0, Rsurf0);

Nsegs=2000;
N_pulses=Nsegs*58;
for k=1:numel(Rough)
    D2_file=sprintf('%s/D2/Rough=%3.2f_Rsurf=%3.2f_BGR=%3.2f.h5', test_data_dir, Rough(k), Rsurf(k), BGR);
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
            [D2(beam), params(beam)]=make_ATL03_data(N_pulses, N_chan(beam), Rough(k), N_per_white(beam)*Rsurf(k), WF, BGR, 40, 0);
            %D2(beam)=index_struct(D2(beam), (D2(beam).h -D2(beam).z0)> -20 & (D2(beam).h-D2(beam).z0) < 30);  % random truncation...          
        end
         
        % assign left-out parameters
        for beam=1:2
            D2(beam).x_RGT=D2(beam).pulse_num*.7;
            D2(beam).y_RGT=zeros(size(D2(beam).x_RGT))+45*beam-90;
            D2(beam).time=D2(beam).pulse_num*1e-4;
            D2(beam).beam=zeros(size(D2(beam).h))+beam;
            D2(beam).segment_number=floor(D2(beam).pulse_num/29)+1;
            % delete parameters that might cause confusion
            x_ps=D2(beam).x_RGT+ll2ps(-82, 116);
            [D2(beam).lat, D2(beam).lon]=ps2ll(x_ps);
        end
        for beam=1:2
            params(beam).N_channels=N_chan(beam);
            params(beam).roughness=Rough(k);
            params(beam).R_surf=Rsurf(k);
        end
           
        write_test_file(D2_file, D2, params);
    end
    if false
        %TEST_CODE% [D2, params]=read_test_file(D2_file);         for beam=1:2; params(beam).WF=WF;end
        %[D2, params]=read_test_file(D2_file);
        for kB=1:2
            D2(kB).t_ph=D2(kB).h/(-c2);
        end
        delta_t=25e-12; %ATBD value
        % make a histogram that's large enough to span the whole range of
        % heights
        t_WF_full=-5e-9:delta_t:5e-9;
        signal_threshold_for_gain_corr=0.01;
        for k0=1:max(D2(2).segment_number)-1
            for kB=1:2
                these=D2(kB).segment_number==k0 | D2(kB).segment_number==k0+1;
                tt=D2(kB).t_ph(these & D2(kB).detected);
                % do an iterative trim on these photons (b/c fpb_corr is
                % tested on trimmed distributions;
                ind=true(size(tt));
                last_ind=false(size(tt));
                it_count=0;
                while any(last_ind~=ind) && it_count<10
                    last_ind=ind;
                    tbar=mean(tt(ind));
                    sigma=max(0.1/c2, iqr(tt(ind))/2);
                    ind=abs(tt-tbar)<3*sigma;
                    it_count=it_count+1;
                end
                
                counts=my_histc(tt(ind)-tbar, t_WF_full);
                [med, centroid, count, N_fpb_corr, sigma_med, sigma_centroid, minGain, gain_full]=...
                    fpb_corr(t_WF_full, counts,  params(kB).N_channels, 58, params(kB).t_dead, signal_threshold_for_gain_corr);
                D3.mean_uncorr(k0, kB)=tbar*(-c2);
                D3.med_uncorr(k0, kB)=median(tt(ind))*(-c2);
                D3.fpb_med_corr(k0, kB)=med;
                D3.fpb_mean_corr(k0, kB)=centroid;
                D3.count(:, k0, kB)=counts;
                D3.count_FPB_corr(:,k0, kB)=N_fpb_corr;
                D3.fpb_sigma_med(k0, kB)=sigma_med;
                D3.fpb_sigma_mean(k0, kB)=sigma_centroid;
                D3.gain(:, k0, kB)=gain_full;
                D3.pulse_0(k0, kB)=min(D2(kB).pulse_num);
                D3.pulse_1(k0, kB)=max(D2(kB).pulse_num);
            end
        end
        D3_file=strrep(D2_file,'D2','D3');
        write_struct_to_h5_file(D3_file, D3);   
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
params=struct('roughness', roughness,   'sigma_x', 7.5, ...
    'NoiseRate', BGR, 'H_window', H_window, 'WF', WF, 'c', 3e8,...
    'N_channels', N_chan, 'N_per_pulse', N_per_pulse,'t_dead', 3.2e-9);
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
        try
            D2(beam).(f2{k})=h5read(filename, ['/photon', beam_name{beam},f2{k}]);
        catch
        end
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

