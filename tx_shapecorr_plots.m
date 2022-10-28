% plots for correct_tx_shape_corr
clear


M(1)=load('tx_shape_corr_test_out_TEP_7_10_2018', 'R','HW' ,'bar_out', 'med_out');
M(2)=load('tx_shape_corr_test_out_WF_7_5_2018', 'R' ,'HW', 'bar_out', 'med_out');
M(1).WF=proc_TEP; M(1).WF.dt=diff(M(1).WF.t(1:2));
load WF_est;WF.dt=diff(WF.t(1:2));
t1=-2e-8:WF.dt:2e-8; WF.p=interp1(WF.t, WF.p, t1); WF.p(~isfinite(WF.p))=0; WF.t=t1;
M(2).WF=WF;
for k=1:2; M(k).WF.p=M(k).WF.p/max(M(k).WF.p); end

YL=[0 1.05];
figure(1); clf; set(gcf, 'units','inches','position', [12 3 6 6], 'defaultaxesfontsize', 12,'color','w');

clear h_ax
h_ax(1:2, 1)=cheek_by_jowl(2,1,[0.15  0.15 0.7 0.7]);
figure(2); clf; set(gcf, 'units','inches','position', [12 3 7 5.2], 'defaultaxesfontsize', 12,'color','w');

h_ax(:, 2:5)=cheek_by_jowl(2,4,[0.10 0.15 0.65 0.7]);
set(h_ax(1:2,1),'yaxislocation','left')
for k=1:2
    axes(h_ax(k,1));
    plot(M(k).WF.t*1e9, M(k).WF.p/max(M(k).WF.p),'k','linewidth', 2);
    hold on;
     
    
    WFc=0;
    els=abs(M(k).WF.t*1.5e8)<1.5;
    for k1=1:10
        els=abs((M(k).WF.t-WFc)*1.5e8)<1.5;
        WFc=sum(M(k).WF.t(els).*M(k).WF.p(els))./sum(M(k).WF.p(els));        
    end
    plot([1 1]*WFc*1e9, YL,'k--','linewidth', 2);
    plot([1 1]*wf_percentile(M(k).WF.t(els)*1e9, M(k).WF.p(els), 0.5), YL,'k:','linewidth', 2);
    Gsamps=ceil(0.25/(M(k).WF.dt*1.5e8))*3;
    WFB=conv(M(k).WF.p, gaussian([-Gsamps:Gsamps]*M(k).WF.dt, 0, 0.25/1.5e8),'same')*M(k).WF.dt;
    plot(M(k).WF.t*1e9, WFB,'r','linewidth', 2);
    WFc=0;
    for k1=1:10
        els=abs((M(k).WF.t-WFc)*1.5e8)<1.5;
        WFc=sum(M(k).WF.t(els).*WFB(els))./sum(WFB(els));        
    end
     
    plot([0 0]+WFc*1e9, YL,'r--','linewidth', 2);
    plot([0 0]+wf_percentile(M(k).WF.t(els), WFB(els), 0.5)*1e9, YL,'r:','linewidth',2);
    set(gca,'xlim', [-6 6]);

    axes(h_ax(k,2))
    imagesc(  M(k).R(1,:), M(k).HW(:,1), M(k).bar_out.med_uncorr);
    hold on;
    plot( M(k).R(1,:), 6*sqrt(M(k).R(1,:).^2+.1^2),  'k--','linewidth', 2);
    axes(h_ax(k,3));
    imagesc( M(k).R(1,:), M(k).HW(:,1), M(k).bar_out.centroid_uncorr);
    hold on;
    plot( M(k).R(1,:), 6*sqrt(M(k).R(1,:).^2+.1^2),  'k--','linewidth', 2);
    
    axes(h_ax(k,4)); 
    imagesc(M(k).R(1,:), M(k).HW(:,1),  M(k).bar_out.med);
    axes(h_ax(k,5));
    imagesc( M(k).R(1,:), M(k).HW(:,1), M(k).bar_out.centroid);
    
end

set(h_ax(1:2,1),'ylim',[0 1.05]);
axes(h_ax(1,1)); 
legend('RX, R=0','mean, R=0','median, R=0','RX, R=0.25','mean, R=0.25','median, R=0.25')
XR=get(gca,'xlim'); YR=get(gca,'ylim');
ht=text(XR(1), YR(2),' a: ATLAS TEP', 'verticalalignment','top','fontsize', 14);

axes(h_ax(2,1)); 
XR=get(gca,'xlim'); YR=get(gca,'ylim');
ht=text(XR(1), YR(2),' b: Prototype', 'verticalalignment','top','fontsize', 14);
set(h_ax(1,1),'xticklabel','');
xlabel(h_ax(2,1),'time, ns');
yl=ylabel(h_ax(1,1),'amplitude (arbitrary units)');
set(yl,'position', get(yl,'position')-[0 .5 0])
 
cax_tep=[-.002 0.012];
cax_wf=[-.012 .12];
% format colored plots:

for k=2:5
    axes(h_ax(1,k)); caxis(cax_tep); set(gca,'ydir','reverse');
    axes(h_ax(2,k)); caxis(cax_wf);  set(gca,'ydir','reverse');
end
set(h_ax,'ydir','normal');

set(h_ax(:, 3:4),'yticklabel','')
set(h_ax(:,5),'yaxislocation','right');
hl=ylabel(h_ax(1,2),'Window width, m');
pos=get(hl,'position');
set(hl,'position', [pos(1), min( M(1).HW(1,:)), pos(3)]);
hl=ylabel(h_ax(1,5),'Window width, m');
pos=get(hl,'position');
set(hl,'position', [pos(1), min( M(1).HW(1,:)), pos(3)]);

hl=xlabel(h_ax(2, 4),'Roughness, m');
pos=get(hl,'position');
set(hl,'position', [min(M(1).R(:,1)), pos(2), pos(3)]);

pos=get(h_ax(1,5),'position');
hb=axes('position', [pos(1)+pos(3)+0.08, pos(2), 0.025, pos(4)]); 
imagesc([0 1], 1000*linspace(cax_tep(1), cax_tep(2), 100), 1000*linspace(cax_tep(1), cax_tep(2), 100)'); set(hb,'xtick'); 
ylabel(hb,{'TEP' 'WF bias, mm'});
set(hb,'yaxislocation','right','ydir','normal');

pos=get(h_ax(2,5),'position');
hb=axes('position', [pos(1)+pos(3)+0.08, pos(2), 0.025, pos(4)]); 
imagesc([0 1], 1000*linspace(cax_wf(1), cax_wf(2), 100), 1000*linspace(cax_tep(1), cax_tep(2), 100)'); set(hb,'xtick'); 
ylabel(hb,{'Prototype' 'WF bias, mm'});
set(hb,'yaxislocation','right','ydir','normal');

title(h_ax(1,2),'median');
title(h_ax(1,3),'mean');
title(h_ax(1,4),{'corrected' 'median'});
title(h_ax(1,5),{'corrected' 'mean'});
colormap(my_rgb_cpt(64)*0.7);
tags={' a',' b',' c',' d';' e',' f',' g',' h'};
XL=get(h_ax(1,2),'xlim');
YL=get(h_ax(1,2),'ylim');
h_tag=[];
for k=2:5
    for k1=1:2;
        axes(h_ax(k1, k)); 
        h_tag(k1,k-1)=text(XL(1), YL(2),tags{k1,k-1}, 'horizontalalignment','left',...
            'verticalalignment','top','fontsize', 14, 'backgroundcolor','w');
    end
    
end

if false
    for k=1:2; pos=get(k,'position'); set(k,'papersize',  pos(3:4)+[0.5 0.5]); end
    print -f1 -dpdf paper_figures/TEP_plus_prototype
    print -f2 -dpdf paper_figures/TX_corr_bias
end
    
    



