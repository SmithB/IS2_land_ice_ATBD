load WF_est;
WF1=WF;
WF.t=(-60:.005:60)*1e-9;
WF.p=interp1(WF1.t, WF1.p, WF.t);
WF.p(~isfinite(WF.p))=0;
WF.t=WF.t-sum(WF.t.*WF.p)/sum(WF.p);
t2h=-1.5e8;

figure(gcf); clf; set(gcf,'defaultaxesfontsize', 12);
subplot(1,2,1)
plot(WF.p/max(WF.p), WF.t*t2h,'k');
hold on;
plot(get(gca,'xlim'), [0 0],'k--');
plot(get(gca,'xlim'), [0 0 ]+ t2h*wf_median(WF.t(abs(WF.t)<1e-8), WF.p(abs(WF.t)<1e-8)),'k:')


P1=conv(WF.p/max(WF.p), gaussian(WF.t,  0, .25/(1.5e8))*diff(WF.t(1:2)) ,'same');
plot(P1, WF.t*t2h,'color', [0.6 0.6 0.6]);

C=sum(WF.t(abs(WF.t)<1e-8).*P1(abs(WF.t)<1e-8))./sum(P1(abs(WF.t)<1e-8));
M=wf_median(WF.t(abs(WF.t)<1e-8), P1(abs(WF.t)<1e-8));
plot(get(gca,'xlim'), [0 0]+C*t2h,'--','color', [0.6 0.6 0.6]);
plot(get(gca,'xlim'), [0 0]+M*t2h,':','color', [0.6 0.6 0.6]);

set(gca,'ylim', [-5.25 5.25]*1e-9*1.5e8);
grid on
legend('Tx','Tx centroid','Tx median','Rx','Rx centroid','Rx median');
xlabel('normalized power'); ylabel('h-Tx centroid, m');

R=0.01:0.01:1.5;
[M1, C1, dM, dCtr, sigma]=deal(NaN(size(R)));
for k=1:length(R)
     P1=conv(WF.p/max(WF.p), gaussian(WF.t,  0, R(k)/(1.5e8))*diff(WF.t(1:2)) ,'same');
     sigma0=diff(wf_percentile(WF.t(:), P1(:), [0.14 0.86]))/2;
     t_win=max(3/1.5e8, 6*sigma0);

     bar_last=-10*sigma0;
     bar=sum(WF.t(abs(WF.t)<t_win/2).*P1(abs(WF.t)<t_win/2))./sum(P1(abs(WF.t)<t_win/2));
     count=0;
     while abs(bar-bar_last) > 0.001*sigma(k) 
         bar_last=bar;
         bar=sum(WF.t(abs(WF.t-bar)<t_win/2).*P1(abs(WF.t-bar)<t_win/2))./sum(P1(abs(WF.t-bar)<t_win/2));
         count=count+1;
     end
     NN(k)=count;
     C1(k)=sum(WF.t(abs(WF.t-bar)<t_win/2).*P1(abs(WF.t-bar)<t_win/2))./sum(P1(abs(WF.t-bar)<t_win/2));
     M1(k)=wf_percentile(WF.t(abs(WF.t-bar)<t_win/2), P1(abs(WF.t-bar)<t_win/2), 0.5);
     
     sigma(k)=diff(wf_percentile(WF.t(abs(WF.t-bar)<t_win/2), P1(abs(WF.t-bar)<t_win/2), [.16 .84]))/2;
     [dM(k), dCtr(k)]=correct_for_TX_shape( WF.t(:), WF.p(:), t_win, sigma(k), 1e6);
end

subplot(1,2,2); cla; hold on;
plot(R, C1*-1.5e8,'b--');
plot(R, C1*-1.5e8+dCtr,'b');
plot(R, M1*-1.5e8,'r--');
plot(R, M1*-1.5e8+dM,  'r');
grid
legend('Centroid','Centroid corr.','Median','Median corr.')
set(findobj(gcf,'type','line'),'linewidth', 2)
xlabel('surface roughness, m'); ylabel('\delta h, m');
 