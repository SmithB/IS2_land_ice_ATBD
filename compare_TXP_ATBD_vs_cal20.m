
csv_file='/home/ben/Dropbox/projects/IS2_ATBD/Pulse_shape_TEP_cal20/6_CAL_20_PRIM_L1_MODE6_2017130_c_M.csv';
[~, out]=unix(['head -71 ', csv_file,' | tail -1']);
out=strsplit(deblank(out),',');
for k=3:length(out)
    WF.p(k-2)=str2num(out{k}); 
    WF.t(k-2)=2.5e-11*(k-2);
end
WF.t=WF.t-sum(WF.t.*WF.p)./sum(WF.p);

H_win=6*[0.1:0.01:3];
for k=1:length(H_win)
    C=0;
    for kk=1:50
        els=abs(WF.t-C)<(H_win(k)/2/1.5e8);
        C=sum(WF.t(els).*WF.p(els))./sum(WF.p(els));
    end   
    TXbar(1,k)=C;
    TXmed(1,k)=wf_median(WF.t(els), WF.p(els)); 
end

load WF_est; 

for k=1:length(H_win)
    C=0;
    for kk=1:50
        els=abs(WF.t-C)<(H_win(k)/2/1.5e8);
        C=sum(WF.t(els).*WF.p(els))./sum(WF.p(els));
    end   
    TXbar(2,k)=C;
    TXmed(2,k)=wf_median(WF.t(els), WF.p(els)); 
end

figure(1); clf; hold on;
set(gca,'colororder', [1 0 0; 0 0 1]);
plot(H_win, TXbar*1.5e8,'o');
plot(H_win, TXmed*1.5e8,'-');


