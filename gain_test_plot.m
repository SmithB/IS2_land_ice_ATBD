function [t_wf, N_det, N_missed]=gain_test_plot(PPC, sigma, NonPara)

if ~exist('PPC','var')
    PPC=1;
end

if ~exist('sigma','var');
    sigma=0.68*1.5;
end

n_ch=16;
n_pulses=57;
t_dead=3.2e-9;
dt=2.5e-10;

n_ph=ceil(n_ch*n_pulses*PPC);
pulse=ceil(rand(n_ph,1)*n_pulses);
ch=ceil(rand(n_ph,1)*n_ch);

t=sigma*randn(n_ph,1);

[~, ind]=sort(pulse*(n_ch+1)+t);
pulse=pulse(ind);
ch=ch(ind);
t=t(ind);
% t_wf=[-3.5e-9:dt:6e-9];
t_wf=[round_to(-2*sigma-t_dead, dt):dt:(3*sigma+t_dead)];
if exist('NonPara','var') && NonPara
    disp('WARNING: USING THE NON-PARALYZABLE CODE')
    [hit_all, det, dead_all]=detect_with_deadtime_nonparalyzable(t, pulse, n_ch, t_dead, n_pulses, t_wf);
    name='nonparalyzable';
else
    [hit_all, det, dead_all]=detect_with_deadtime(t, pulse, n_ch, t_dead, n_pulses, t_wf);
    name='   paralyzable';
end
n_all=histcounts(t, t_wf);
n_detected=histcounts(t(hit_all), t_wf);

dt_fine=2.5e-12;
t_wf_fine=t_wf(1):dt_fine:t_wf(end);
n_detected_fine=my_histc(t(hit_all), t_wf_fine);
[med, centroid, count, N_fpb_corr, sigma_med, sigma_centroid, minGain, gain_fine]=...
    fpb_corr(t_wf_fine, n_detected_fine,  n_ch, n_pulses, t_dead, 1e-3, dt_fine);
gain_full=interp1(t_wf_fine+dt_fine/2, gain_fine, t_wf+dt/2);


figure; clf; 
set(gcf, 'units','inches','position', [12 3 6 6], 'defaultaxesfontsize', 12,'color','w','papersize', [6.5 6.5]);
hold on; 

hb1=bar((t_wf(1:end-1)+dt/2)*1e9, n_all/max(n_all),1,'r'); set(hb1,'facecolor', [ 1 0 0 ]*.8+.1,'edgecolor','none')
hb2=bar((t_wf(1:end-1)+dt/2)*1e9, n_detected/max(n_all),1,'b'); set(hb2,'facecolor', [ 0 0 1 ]*.8+.1,'edgecolor','none')

hl1=plot((t_wf(1:end-1)-dt/2)*1e9, 1-mean(dead_all(1:end-1,:),2),'k','linewidth', 3);
hl2=plot(t_wf*1e9, gain_full,'color', [1 1 1]*.7,'linewidth',3);

hb3=bar((t_wf(1:end-1)+dt/2)*1e9, n_detected./gain_full(1:end-1)/max(n_all)); set(hb3,'facealpha', 0.5,'facecolor', [1 1 1]*.5)
legend([hb1, hb2, hl1, hl2, hb3],{'true count','detected count','true gain','est gain','est. count'});

xlabel('time, ns');
ylabel('normalized count, gain')
if false
    print -f1 -dpdf paper_figures/gain_est
end
fprintf(1,'%s, mean_det = %3.2f, est_mean_det=%3.2f, med_z=%3.2f mean_z=%3.2f\n', name, mean(hit_all), sum(N_fpb_corr)/length(hit_all),  median(t(hit_all))*1.5e8, mean(t(hit_all))*1.5e8)
fprintf(1,'\tNPP_true=%2.1f, NPP_det=%2.1f, NPP_corr=%2.1f\n', n_ph/n_pulses, sum(hit_all)/n_pulses, sum(N_fpb_corr)/n_pulses);




%----------------------------------------------------
function [hit_all, det, dead_all] = detect_with_deadtime_nonparalyzable(t, pulse,  N_det, t_dead, N_pulses, t_WF)
% NOTE: t and pulse must be sorted by pulse, then by time.
DS=ceil(rand(size(t))*N_det)+1i*pulse;
u_DS=unique(DS);

if exist('t_WF','var')
    dead_all=zeros(length(t_WF), N_pulses*N_det);
    dt=t_WF(2)-t_WF(1);
    N_dead_bins=ceil(t_dead/dt);
    calc_gain=true;
else
    calc_gain=false;
end

hit_all=false(size(t));
for k=1:length(u_DS)
    these=DS==u_DS(k);
    if ~any(these) 
        continue; 
    end
    % ti are the photons for this detector
    ti=t(these);
    hit=false(size(ti));
    missed=false(size(ti));
    next_first=1;
    hit(next_first)=true;
    
    if calc_gain
        dead=zeros(size(t_WF));
        dead(t_WF>ti(1) & t_WF < ti(1)+t_dead) = true;
    end
    %disp(k)
    while ~(hit(end) || missed(end))
        first=next_first;
        hit(first)=true;
        in_dt=ti>ti(first) & ti < ti(first)+t_dead;
        if any(in_dt)
            % extend the dead time if there are photon hits while the counter is dead
            count=0;
             
            missed(in_dt) = true;
            last=find(in_dt, 1, 'last');
            if calc_gain
                dead(t_WF>ti(first) & t_WF < ti(first)+t_dead) = true;
            end
            next_first=last+1;
        else
            if calc_gain
                dead(t_WF>ti(first) & t_WF < ti(first)+t_dead) = true;
            end
            next_first=first+1;
        end
    end
    if calc_gain
        dead_all(:, real(u_DS(k))+imag(u_DS(k))*N_det)=dead(:);
    end
    hit_all(these)=hit;
end
det=real(DS);










%----------------------------------------------------
function [hit_all, det, dead_all] = detect_with_deadtime(t, pulse,  N_det, t_dead, N_pulses, t_WF)
% NOTE: t and pulse must be sorted by pulse, then by time.
DS=ceil(rand(size(t))*N_det)+1i*pulse;
u_DS=unique(DS);

if exist('t_WF','var')
    dead_all=zeros(length(t_WF), N_pulses*N_det);
    dt=t_WF(2)-t_WF(1);
    N_dead_bins=ceil(t_dead/dt);
    calc_gain=true;
else
    calc_gain=false;
end

hit_all=false(size(t));
for k=1:length(u_DS)
    these=DS==u_DS(k);
    if ~any(these) 
        continue; 
    end
    % ti are the photons for this detector
    ti=t(these);
    hit=false(size(ti));
    missed=false(size(ti));
    next_first=1;
    hit(next_first)=true;
    
    if calc_gain
        dead=zeros(size(t_WF));
        dead(t_WF>ti(1) & t_WF < ti(1)+t_dead) = true;
    end
    %disp(k)
    while ~(hit(end) || missed(end))
        first=next_first;
        hit(first)=true;
        in_dt=ti>ti(first) & ti < ti(first)+t_dead;
        if any(in_dt)
            % extend the dead time if there are photon hits while the counter is dead
            count=0;
            while sum(in_dt) > count
                count=sum(in_dt);
                in_dt=ti > ti(first) & ti < max(ti(in_dt))+t_dead;
            end
            missed(in_dt) = true;
            last=find(in_dt, 1, 'last');
            if calc_gain
                dead(t_WF>ti(first) & t_WF < ti(last)+t_dead) = true;
            end
            next_first=last+1;
        else
            if calc_gain
                dead(t_WF>ti(first) & t_WF < ti(first)+t_dead) = true;
            end
            next_first=first+1;
        end
    end
    if calc_gain
        dead_all(:, real(u_DS(k))+imag(u_DS(k))*N_det)=dead(:);
    end
    hit_all(these)=hit;
end
det=real(DS);
