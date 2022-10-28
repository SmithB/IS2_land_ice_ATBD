function plot_nick_profile(s, D)

s=str2num(s);
h_ax(1)=axes('position', [0.15 0.4 0.7 0.6]);
hold on;
plot(D(s(1), s(2)).x_atc, D(s(1), s(2)).DEM,'k.','markersize', 1);
good=D(s(1), s(2)).atl06_quality_summary==0;
for kB=1:2
     
    plot(D(s(1), s(2)).x_atc(good(:, kB), kB), D(s(1), s(2)).h_li(good(:, kB), kB),'.');
end
h_ax(2)=axes('position', [0.15 0.15 0.7 0.25]);
plot(D(s(1), s(2)).x_atc, conv2(double(good), ones(50,1)/50,'same'),'.');

linkaxes(h_ax,'x');
