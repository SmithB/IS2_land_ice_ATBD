function WF=read_cal20_wf

csv_file='/home/ben/Dropbox/projects/IS2_ATBD/Pulse_shape_TEP_cal20/6_CAL_20_PRIM_L1_MODE6_2017130_c_M.csv';
[~, out]=unix(['head -71 ', csv_file,' | tail -1']);
out=strsplit(deblank(out),',');
for k=3:length(out)
    WF.p(k-2)=str2num(out{k}); 
    WF.t(k-2)=2.5e-11*(k-2);
end
WF.t=WF.t-sum(WF.t.*WF.p)./sum(WF.p);

