% set a breakpoint at test_one_SNR_significance_val.m, line 40, then run
% this.

figure(10); clf; set(gcf,'defaultaxesfontsize', 12,'units','inches','position', [1 1 5.5 4.5]);

clear FF;
these=selected_PE;

plot(D2sub.x_RGT, D2sub.h,'ko','markersize', 4,'markerfacecolor','k'); hold on;

plot(D2sub.x_RGT(these), D2sub.h(these),'ro');
xlabel('along-track, m'); ylabel('height, m');
N=sum(these); 
HW=diff(range(D2sub.h(these)));
H0=mean(range(D2sub.h(these)));
dhdx=0;
N_noise=1e7/1.5e8*HW; SNR=(N-N_noise)/N_noise;
title(sprintf('N=%d, W=%2.1f, SNR=%2.2f', sum(selected_PE), HW, SNR));
set(gca,'xlim', x_seg_ctr+[-20 20],'ylim', [-10 10]);
FF(1)=getframe(gcf);
for k=1:length(LOG.iterations)
    %figure(k+10);
    clf;
    plot(D2sub.x_RGT, D2sub.h,'ko','markersize', 6,'markerfacecolor','k'); hold on;
    plot(D2sub.x_RGT(these), D2sub.h(these),'ro','linewidth', 2);
    
    %plot_segs(x_seg_ctr, H0, dhdx,40,{'color', [.7 .7 .7],'linewidth', 3});
    plot_segs(x_seg_ctr, H0-HW/2, dhdx,40,{'--','color', [.7 .7 .7],'linewidth', 1});
    plot_segs(x_seg_ctr, H0+HW/2, dhdx,40,{'--','color', [.7 .7 .7],'linewidth', 1});
    
    set(gca,'xlim', x_seg_ctr+[-20 20],'ylim', [-10 10]);
    xlabel('along-track, m'); ylabel('height, m');

    FF(end+1)=getframe(gcf);

    
    these=LOG.iterations(k).els;
    N=sum(these); HW=LOG.iterations(k).W_win;
    H0=LOG.iterations(k).h_ctr;
    dhdx= LOG.iterations(k).dhdx;
    plot_segs(x_seg_ctr, H0, dhdx,40,{'r','linewidth', 3});
    %plot_segs(x_seg_ctr, H0-HW/2, dhdx,40,{'r--','linewidth', 1});
    %plot_segs(x_seg_ctr, H0+HW/2, dhdx,40,{'r--','linewidth', 1});
        
    N_noise=1e7/1.5e8*HW; SNR=(N-N_noise)/N_noise;
    title(sprintf('N=%d, W=%2.1f, SNR=%2.2f', sum(these), HW, SNR));
    set(gca,'xlim', x_seg_ctr+[-20 20],'ylim', [-10 10]);
    FF(end+1)=getframe(gcf);
end

for k=1:5; FF(end+1)=getframe(gcf); end

VV=VideoWriter('SigFinder.avi'); VV.FrameRate=1; 
open(VV);
for k=1:length(FF); writeVideo(VV, FF(k)); end;
close(VV);

