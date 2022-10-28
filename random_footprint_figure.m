[xxs,yys]=meshgrid((0:.7:40)-20, 45+zeros(12,1));
xgs=xxs+randn(size(xxs))*7.5+1i*(yys+7.5*randn(size(yys)));
zs=pi/180*(imag(xgs)-45); 
zs=zs+randn(size(zs))*0.68*.15;

[xxw,yyw]=meshgrid((0:.7:40)-20, -45+zeros(3,1));
xgw=xxw+randn(size(xxw))*7.5+1i*(yyw+7.5*randn(size(yyw)));
zw=pi/180*(imag(xgw)+45); 
zw=zw+randn(size(zw))*0.68*.15;


figure('units','inches','position', [1 1 6 7.5],'defaultaxesfontsize', 12);
axes('position', [0.05 0.1 0.7 0.8]); set(gca,'xlim', [-40 40],'ylim', [-75 75],'dataaspectratio', [1 1 1 ],'fontsize', 12)
hold on;
plot([-40 40], [45 45],'k--'); plot([-40 40], [-45 -45],'k--');
plot(xgs(:,29),'g.','markersize', 15); plot(xgw(:,29),'g.','markersize', 15);
plot(xxs(1, 29), yys(1,29),'ko');
plot(xxw(1, 29), yyw(1,29),'ko');

xlabel('along-track, m'); ylabel('across-track, m')
h1=axes('position', [0.7, 0.60 0.2, 0.25],'fontsize', 12,'yaxislocation','right');
histogram(zs(:,29), -1:.05:1); xlabel('\delta z, m'); ylabel('count'); title('N=12');
set(gca,'yaxislocation','right');

h2=axes('position', [0.7, 0.15 0.2, 0.25],'fontsize', 12); 
histogram(zw(:,29), -1:.05:1); xlabel('\delta z, m'); ylabel('count'); title('N=3');
set(gca,'yaxislocation','right');


figure('units','inches','position', [1 1 6 7.5 ],'defaultaxesfontsize', 12);
axes('position', [0.05 0.1 0.7 0.8]); set(gca,'xlim', [-40 40],'ylim', [-75 75],'dataaspectratio', [1 1 1 ],'fontsize', 12)
hold on;
plot([-40 40], [45 45],'k--'); plot([-40 40], [-45 -45],'k--');
plot(xgs(:),'g.','markersize', 12); plot(xgw(:),'g.','markersize', 12);
plot(xgs(:,29),'r.','markersize', 15); plot(xgw(:,29),'r.','markersize', 15);
plot(xxs(1, :), yys(1,29),'ko');
plot(xxw(1, :), yyw(1,29),'ko');

xlabel('along-track, m'); ylabel('across-track, m')
h1=axes('position', [0.7, 0.60 0.2, 0.25],'fontsize', 12,'yaxislocation','right');
histogram(zs(:), -1:.05:1); xlabel('\delta z, m'); ylabel('count'); title('N=696');
set(gca,'yaxislocation','right');

h2=axes('position', [0.7, 0.15 0.2, 0.25],'fontsize', 12,'yaxislocation','right');
histogram(zw(:), -1:.05:1); xlabel('\delta z, m'); ylabel('count'); title('N=174');
set(gca,'yaxislocation','right');


 