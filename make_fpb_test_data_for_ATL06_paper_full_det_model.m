function varargout=make_fpb_test_data_for_ATL06_paper(varargin)

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

Rough0=([ 0, 2 4, 6 8  10 12  ])*sigma0;
Rsurf0=[0.125:0.125:1.5];
BGR=0;%1e5;%1e7;
test_data_dir=['/home/ben/temp/deadtime_correction_data/tdead_dig=3.2ns.t_dead_analog=1ns'];
out_dir=test_data_dir; 

N_chan=[4 16];
N_per_white=N_chan*3/4;

if ~exist(test_data_dir,'dir')
    mkdir(test_data_dir);
end
if ~exist([test_data_dir,'/D2'],'dir')
    mkdir([test_data_dir,'/D2']);
end

if ~exist([test_data_dir,'/D3'],'dir')
    mkdir([test_data_dir,'/D3']);
end

[Rough,Rsurf]=meshgrid(Rough0, Rsurf0);

Nsegs=50000;
H_win=40;
N_pulses=Nsegs*58;
parfor k=1:numel(Rough)
    make_one_test_file(test_data_dir, Rough, Rsurf, k, BGR, N_pulses, N_chan, H_win, N_per_white, WF)
end

%-------------------------------------------------------------------------------------------------
function make_one_test_file(test_data_dir, Rough, Rsurf, k, BGR, N_pulses, N_chan, H_win, N_per_white, WF)
c2=1.5e8;
D2_file=sprintf('%s/D2/Rough=%3.2f_Rsurf=%3.2f_BGR=%3.2f.h5', test_data_dir, Rough(k), Rsurf(k), BGR);
fprintf(1,'D2_file=%s\n', D2_file);
if exist(D2_file,'file')
    [D2, params]=read_test_file(D2_file);
    for beam=1:2; params(beam).WF=WF;end
else
    clear D2 params;
    for beam=1:2
        [D2(beam), params(beam)]=make_ATL03_data(N_pulses, N_chan(beam), Rough(k), N_per_white(beam)*Rsurf(k), WF, BGR, H_win, 0);
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

for kB=1:2
    D2(kB).t_ph=D2(kB).h/(-c2);
end
delta_t=25e-12; %ATBD value
% make a histogram that's large enough to span the whole range of
% heights
t_WF_full=-5e-8:delta_t:5e-8;
signal_threshold_for_gain_corr=0.01;
%[D3.count_uncorrected, D3.count_FPB_corr, D3.gain]=deal(NaN(length(t_WF_full), max(D2(2).segment_number)-1, 2));
[D3.N_true, D3.N_uncorr, D3.mean_uncorr, D3.med_uncorr, D3.fpb_med_corr, D3.fpb_mean_corr, D3.N_corrected, D3.fpb_sigma_med, D3.fpb_sigma_mean, D3.pulse_0, D3.pulse_1]=...
    deal(NaN(max(D2(2).segment_number)-1, 2));

for kB=1:2
    els{kB}=bin_by(D2(kB).segment_number, 1:max(D2(kB).segment_number));
end

tic;
for k0=1:max(D2(2).segment_number)-1   
    for kB=1:2
        try
            %these=D2(kB).segment_number==k0 | D2(kB).segment_number==k0+1;
            these=[els{kB}{k0}; els{kB}{k0+1}];
            tt_all=D2(kB).t_ph(these);
            tt=tt_all(D2(kB).detected(these)==1);
            N_pulses_est=max(57,diff(range(D2(kB).pulse_num(these)))+1);
            % do an iterative trim on these photons (b/c fpb_corr is
            % tested on trimmed distributions;
            ind=true(size(tt));
            last_ind=false(size(tt));
            it_count=0;
            sigma=H_win/6/c2;
            sigma_last=sigma;
            while any(last_ind~=ind) && it_count<50
                last_ind=ind;
                tbar=mean(tt(ind));
                sigma_last=sigma;
                if BGR==0
                    sigma=max(0.5/c2, iqr(tt(ind))/2);
                else
                    [sigma_h]=robust_peak_width_CDF(tt*-c2, 6*sigma*BGR*N_pulses_est, [-3 3]*sigma);
                    sigma=max([0.5/c2, sigma_h/c2, 0.8*sigma_last]);
                end
                ind=abs(tt-tbar)<3*sigma;
                it_count=it_count+1;
            end

            D3.N_true(k0, kB)=sum(abs(tt_all-tbar)<3*sigma);
            D3.N_uncorr(k0,kB)=sum(ind);

            uncorrected_counts=my_histc(tt(ind)-tbar, t_WF_full);
            [med, centroid, count, N_fpb_corr, sigma_med, sigma_centroid, minGain, gain_full]=...
                fpb_corr_nonparalyzable(t_WF_full, uncorrected_counts,  params(kB).N_channels, N_pulses_est, params(kB).t_dead_digital, delta_t);
            D3.mean_uncorr(k0, kB)=mean(tt(ind))*(-c2);
            D3.med_uncorr(k0, kB)=median(tt(ind))*(-c2);
            D3.fpb_med_corr(k0, kB)=med;
            D3.fpb_mean_corr(k0, kB)=centroid;
            %D3.count_uncorrected(:, k0, kB)=uncorrected_counts;
            D3.N_corrected(k0, kB)=sum(N_fpb_corr);
            %D3.count_FPB_corr(:,k0, kB)=N_fpb_corr;
            D3.fpb_sigma_med(k0, kB)=sigma_med;
            D3.fpb_sigma_mean(k0, kB)=sigma_centroid;
            %D3.gain(:, k0, kB)=gain_full;
            D3.pulse_0(k0, kB)=min(D2(kB).pulse_num(these));
            D3.pulse_1(k0, kB)=max(D2(kB).pulse_num(these));
        catch 
            disp(lasterr)
            fprintf(1, 'problem with %s at segment %d beam %d\n' , D2_file, k0, kB);
        end
    end
    if mod(k0, 5000)==0
         t_now=toc;
         fprintf(1, '%s: k0=%d,  %4.0f s elapsed out of %6.0f s expected\n', D2_file, k0, t_now, t_now/k0*max(D2(1).segment_number));
    end
end
D3_file=strrep(D2_file,'D2','D3');
write_struct_to_h5_file(D3_file, D3);

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
    'N_channels', N_chan, 'N_per_pulse', N_per_pulse,'t_dead_digital', 3.2e-9, 't_dead_analog', 1e-9);
D2=det_sim_full_det_model(DEM, x0, params, ATM_xmit);


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
fp={'H_window','N_channels','c','N_per_pulse','t_dead_analog', 't_dead_digital', 'sigma_x','roughness','R_surf'};

%-----------------------------------------------
function write_test_file(filename, D2,  params)

beam_name={'/weak/','/strong/'};
[f2, fp]=define_fields;
for k=1:length(f2)
    for beam=1:2
        if isfield(D2(beam),f2{k})
            write_h5_field(filename,  D2(beam).(f2{k}), ['/photon',beam_name{beam},f2{k}]);
        else
            %warning('%s is not a field of D2', f2{k});
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
        try
            params(beam).(fp{k})=h5read(filename, ['/params', beam_name{beam},fp{k}]);
        catch
        end
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

%------------------------------------------------------
function [Rough, Rsurf, U, C]=report_errors(thedir)
%[Rough, Rsurf, U, C]=make_fpb_test_data('report_errors',thedir)

fields={'fpb_mean_corr','fpb_med_corr','mean_uncorr','med_uncorr' ,'N_uncorr','N_corrected'};

[~, out]=unix(['ls ', thedir]);
out=strsplit(deblank(out));
for kf=1:length(out)
    temp=regexp(out{kf}, 'Rough=(.*)_Rsurf=(.*)_BGR','tokens');
    Rough(kf)=str2num(temp{1}{1});
    Rsurf(kf)=str2num(temp{1}{2});
    thefile=[thedir, '/',out{kf}];
    for k=1:length(fields);
        D3.(fields{k})=h5read(thefile, ['/', fields{k}]);
    end
    for kB=1:2
        U(kB).mean(kf)=mean(D3.mean_uncorr(:,kB));
        U(kB).med(kf)=mean(D3.med_uncorr(:, kB));
        U(kB).N(kf)=mean(D3.N_uncorr(:, kB));
        U(kB).sigma(kf)=std(D3.med_uncorr(:, kB));
        
        C(kB).mean(kf)=mean(D3.mean_uncorr(:,kB)+D3.fpb_mean_corr(:,kB));
        C(kB).med(kf)=mean(D3.mean_uncorr(:,kB)+D3.fpb_med_corr(:,kB));
        C(kB).N(kf)=mean(D3.N_corrected(:, kB));
        C(kB).sigma(kf)=std( D3.mean_uncorr(:,kB)+D3.fpb_med_corr(:,kB));
    end
end
[~, ~, row]=unique(Rough);
[~, ~, col]=unique(Rsurf);

F=fieldnames(U);
temp=zeros(max(row), max(col));
for kF=1:length(F)
    for kB=1:2
        temp(sub2ind(size(temp), row, col))=U(kB).(F{kF});
        U(kB).(F{kF})=temp;
    end
end
F=fieldnames(C);
for kF=1:length(F)
    for kB=1:2
        temp(sub2ind(size(temp), row, col))=C(kB).(F{kF});
        C(kB).(F{kF})=temp;
    end
end

temp(sub2ind(size(temp), row, col))=Rough;
Rough=temp;
temp(sub2ind(size(temp), row, col))=Rsurf;
Rsurf=temp;

%--------------------------------------
function make_figure
thedir=['/home/ben/temp/deadtime_correction_data/tdead_dig=3.2ns.t_dead_analog=1ns/D3/'];
%thedir='/Volumes/insar7/ben/deadtime_correction_data/tdead=3.2_0MHz/D3/';
[Rough, Rsurf, U, C]=make_fpb_test_data_for_ATL06_paper('report_errors', thedir);
figure(51); clf;
 
cmap=my_rgb_cpt(128)*.8;
cmap=cmap([1:2:64, 65:128], :);
colormap(cmap)
set(gcf, 'units','inches','position', [12 3 7 5.2], 'papersize', [7.5 5.5], 'defaultaxesfontsize', 14,'color','w','defaultaxesydir','normal');
hax=cheek_by_jowl(2,3,[.10 0.1 0.8 0.8]);
axes(hax(1,1)); 
imagesc(Rsurf(1,:), Rough(:,1), U(2).mean*1000); caxis([-2 44]);
title('mean');
axes(hax(2,1)); 
imagesc(Rsurf(1,:), Rough(:,1), C(2).mean*1000); caxis([-2 44]);
hb(1)=colorbar('north'); set(hb(1),'color', 'w'); ylabel(hb(1),'bias, mm');

axes(hax(1,2)); 
imagesc(Rsurf(1,:), Rough(:,1), U(2).med*1000); caxis([-2 44]);
title('median');
axes(hax(2,2)); 
imagesc(Rsurf(1,:), Rough(:,1), C(2).med*1000); caxis([-2 44]);
hb(2)=colorbar('north'); set(hb(2),'color', 'w'); ylabel(hb(2),'bias, mm');

axes(hax(1,3));
imagesc(Rsurf(1,:), Rough(:,1), 100*U(2).N/(57*12)./(repmat(.125:.125:1.5, size(Rough,1), 1))-100); caxis([0.65 1.05]*100-100);
colormap(gca,'gray');
title('reflectance');
axes(hax(2,3));
imagesc(Rsurf(1,:), Rough(:,1), 100*C(2).N/(57*12)./(repmat(.125:.125:1.5, size(Rough,1), 1))-100); caxis([0.65 1.05]*100-100);
colormap(gca,'gray');
hb(3)=colorbar('north');
ylabel(hb(3),'reflectance error, %');
set(hax(:,2),'ytick',[]);
set(hax(:,3),'yaxislocation','right');
set(hax,'ydir','normal');
xlabel(hax(2,2),'surface reflectance');
yl=ylabel(hax(2,1),'roughness, m'); pos=get(yl,'position'); set(yl,'position', [pos(1), 1.2, 1]);
yl=ylabel(hax(2,3),'roughness, m'); pos=get(yl,'position'); set(yl,'position', [pos(1), 1.2, 1]);

ht=text(0.1, 1.15,'a','backgroundcolor','w','parent', hax(1,1),'fontsize', 16);
ht=text(0.1, 1.15,'b','backgroundcolor','w','parent', hax(1,2),'fontsize', 16);
ht=text(0.1, 1.15,'c','backgroundcolor','w','parent', hax(1,3),'fontsize', 16);

ht=text(0.1, 1.15,'d','backgroundcolor','w','parent', hax(2,1),'fontsize', 16);
ht=text(0.1, 1.15,'e','backgroundcolor','w','parent', hax(2,2),'fontsize', 16);
ht=text(0.1, 1.15,'f','backgroundcolor','w','parent', hax(2,3),'fontsize', 16);

for k=1:3
    pos=get(hb(k),'position'); 
    set(hb(k),'position', [pos(1)+pos(3)*.25, pos(2), pos(3)*.75, pos(4)]);
end
set(hax(1,:),'xticklabel','');

