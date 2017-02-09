function t=random_wf_sim(WF_t, WF_p, n_photons, BG0, BG1);


dt=(WF_t(end)-WF_t(1))/(length(WF_t)-1);

if exist('BG0','var');
    WF_p=max(0, WF_p-BG0*dt);
end

if exist('BG1','var');
    WF_p=max(0, WF_p+BG1*dt);
end

C_t=[WF_t(1)-dt/2; WF_t(:)+dt/2];
C=[0; cumsum(WF_p(:).*dt)];
C=C/C(end);
if size(n_photons)==1
    r=rand(n_photons,1);
else
    r=rand(n_photons);
end

% generate a sorted version of C
[C1, iC1]=unique(C,'first');
[~, iC2]=unique(C,'last');
C_ti=(C_t(iC1)+C_t(iC2))/2;

t=interp1(C1, C_ti, r,'*linear');


