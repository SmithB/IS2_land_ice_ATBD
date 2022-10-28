function [t_wf, N_det, N_missed]=gain_test_plot_2dtmodels(PPC_vals);
figure; clf; 
set(gcf, 'units','inches','position', [12 3 6 6], 'defaultaxesfontsize', 12,'color','w','papersize', [6.5 6.5]);

if ~exist('PPC_vals','var')
    PPC_vals=[1 2 4 8 16];
end

hax=cheek_by_jowl(1, length(PPC_vals), [0.1 0.1 0.8 0.8]);

for k_ppc=1:length(PPC_vals)
    PPC=PPC_vals(k_ppc);
    n_ch=16;
    n_pulses=5700;
    t_dead_analog=1e-9;
    t_dead_digital=3.2e-9;
    
    % histogram bin size
    dt=2.5e-10;
    
    n_ph=ceil(n_ch*n_pulses*PPC);
    pulse=ceil(rand(n_ph,1)*n_pulses);
    ch=ceil(rand(n_ph,1)*n_ch);
    
    load WF_SuperTep;
    %t=sigma*randn(n_ph,1);
    t=random_wf_sim(WF.t, WF.p, n_ph);
    sigma=0.68e-9;
    
    [~, ind]=sort(pulse*(n_ch+1)+t);
    pulse=pulse(ind);
    ch=ch(ind);
    t=t(ind);
    t_wf=[-5e-9:dt:3e-8];
    dt_params=struct('t', t_wf,'dt', dt,'n_pulses', n_pulses,'pulse_0', 1);
    [hit_all, det, dead_all] = deadtime2(t, pulse,  n_ch, t_dead_analog, t_dead_digital,  dt_params);
    %[hit_all, det, dead_all]=detect_with_deadtime_nonparalyzable(t, pulse, n_ch, t_dead, n_pulses, t_wf);
    
    n_all=histcounts(t, t_wf);
    n_detected=histcounts(t(hit_all), t_wf);
 
    axes(hax(k_ppc));
    hb1=barh((t_wf(1:end-1)+dt/2)*-1.5e8, n_detected/max(n_detected),1,'k'); set(hb1,'facecolor', [ 1 0 0 ]*.8+.1,'edgecolor','none')
    set(gca,'xlim', [0 0.05]);
    title(num2str(PPC));
end


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
function [hit_all, det, dead_all] = detect_with_deadtime_paralyzable(t, pulse,  N_det, t_dead, N_pulses, t_WF)
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
