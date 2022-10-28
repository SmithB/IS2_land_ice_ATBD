thedir='/home/ben/Dropbox/projects/IS2_ATBD/Combined_corr_test_data//September_18_2018_new_WF/D3/';
load([thedir,'/Rough=1.20e-01_Rsurf=6.25e-02_BGR=1.00e+07_LOG.mat'])
%load Rough=9.60e-01_Rsurf=1.25e-01_BGR=1.00e+07_LOG.mat
dhdx_true=.0125;


for k=1:40
    fprintf(1,'%d\t', k)
    for kB=1:2;
        fprintf(1, '%d\t%3.2f \t\t', LOG(k,kB).D3.exit_iteration, LOG(k,kB).D3.h_mean-dhdx_true*LOG(k,kB).D3.x_RGT)
    end
    fprintf(1,'\n');
end

seg0=[3 14 3];
beam0=[2 1 1];



figure(1); clf;
set(gcf,'defaultaxesfontsize', 11,'units','inches','paperpositionmode','auto','position', [1 1 8 4],'papersize', [8.5 4.5],'color','w','inverthardcopy','off');
hax=cheek_by_jowl(1,length(beam0), [0.1 0.15 0.8 0.7], [0.05 0.01]);
titles={'a','b','c'};

for kk=1:length(seg0)
    %for kSeg=1:length(seg0)
    %    seg=seg0(kSeg);
    %     for kB=1:2
    %        beam=beam0(kB);
    seg=seg0(kk);
    beam=beam0(kk);
    LOGi=LOG(seg, beam);
    
    %axes(hax(kB));
    axes(hax(kk));
    cla; hold on;
    x0=LOGi.D3.x_RGT;
    h0=x0*dhdx_true;
    
    colors=jet(14);
    h_seg=[];
    legT={};
    for ki=1:length(LOGi.LS_fit.iterations)
        temp=LOGi.LS_fit.iterations(ki);
        h_seg(ki)=plot([-20 20], temp.h_ctr+temp.dhdx*[-20 20]-h0,'color', colors(ki,:),'linewidth', 2.5);
        plot([-20 20], temp.h_ctr+temp.dhdx*[-20 20]-temp.W_win/2-h0,'--','color', colors(ki,:),'linewidth', 2.5);
        plot([-20 20], temp.h_ctr+temp.dhdx*[-20 20]+temp.W_win/2-h0,'--','color', colors(ki,:),'linewidth', 2.5);
        legT{ki}=sprintf('iteration %d', ki);
    end
    if length(h_seg) > 1
        h_seg=h_seg([1 end]);
        legT=legT([1 end]);
    end
    h_seg(end+1)=plot(LOGi.D2.x_RGT-x0, LOGi.D2.h-h0,'ko', 'linewidth', 1.5); hold on;
    legT{end+1}='all ph.';
    h_seg(end+1)=plot(LOGi.D2.x_RGT(temp.els)-x0, LOGi.D2.h(temp.els)-h0,'ko', 'linewidth', 2,'markerfacecolor','k'); hold on;
    legT{end+1}='sel. ph.';
    h_seg(end+1)=plot([-20 20], dhdx_true*x0+[-20 20]*dhdx_true-h0,'k','linewidth', 2);
    legT{end+1}='surface';
    if kk==1
        legend(h_seg([1  end]), legT([1  end]))
    else
        legend(h_seg([1 2 end]), legT([1 2 end]))
    end
    title(sprintf('$\\mathsf{%s: N=%d, H_{win}=%1.1fm}$,\n $\\mathsf{SNR=%0.1f, SNR_{sig.}=%2.2f}$', titles{kk}, sum(LOGi.LS_fit.iterations(end).els), LOGi.D3.w_surface_window_final, LOGi.D3.SNR, LOGi.D3.SNR_significance),'interpreter','latex');
end
set(hax,'xlim', [-19 19],'ylim', [-19 19]);

set(hax(2),'yticklabel','');

xlabel(hax(2),'x_{atc} (m)');
ylabel(hax(1),'height (m)');
ylabel(hax(3),'height (m)');

%print -dpdf ~/Dropbox/ATL06_paper/paper_figures/SigRefinementExamples.pdf
