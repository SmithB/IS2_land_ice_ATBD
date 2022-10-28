function ATL11_plot(D3b)
clf
hax1=cheek_by_jowl(2,1, [0.1 0.4 0.8 0.5]);
set(hax1,'yaxislocation','left');
hax2=[subplot(3,2,5), subplot(3,2,6)];
%h(:,2)=cheek_by_jowl(3,1, [0.62 0.075 0.34 0.85]);

axes(hax1(1));
colors=jet(12)*0.8;
for k=1:size(D3b.corrected_h.mean_pass_time,2);
    good=D3b.pass_quality.summary(:,k)==0;
    plot(D3b.reference_point.x_ATC, D3b.corrected_h.pass_h_shapecorr(:,k),'.','color', colors(k,:),'markersize', 1); hold on;
    plot(D3b.reference_point.x_ATC(good), D3b.corrected_h.pass_h_shapecorr(good,k),'.','color', colors(k,:),'markersize', 10); hold on;

    %M=[NaN -1  1 NaN]'*D3b.corrected_h.pass_h_shapecorr_sigma(:,k)' + [NaN 1 1 NaN]'*D3b.corrected_h.pass_h_shapecorr(:,k)';
    %plot(ones(4,1)*D3b.corrected_h.ref_pt_X_RPT', M,'color', colors(k,:));
end

%axes(h(2));
%good=D3b.pass_quality.summary==0;
% hL=plot(D3b.reference_point.x_ATC, D3b.ref_surf.N_cycle_avail,'k', D3b.reference_point.x_ATC, sum(good,2),'r');
% set(hL,'linewidth', 2);
% legend(hL,'available','passed');
% ylabel('rep count');
 
Nrep= size(D3b.corrected_h.pass_h_shapecorr,2);
Npt= size(D3b.corrected_h.pass_h_shapecorr,1);
W=D3b.corrected_h.pass_h_shapecorr_sigma.^2+D3b.corrected_h.pass_h_shapecorr_sigma_systematic.^2;
W(D3b.pass_quality.summary>0)=0;
W(W~=0)=1./W(W~=0); W(~isfinite(W))=0;
sumW=sum(W,2);
W(sumW>0,:)=W(sumW>0,:)./repmat(sumW(sumW>0), [1, Nrep]);
D3b.ref_surf.h_corr_mean=NaN(Npt,1);
if ~isfield(D3b.non_product,'delta_h_true');
    D3b.corrected_h.delta_h_true=zeros(size(D3b.corrected_h.pass_h_shapecorr));
end
temp=D3b.corrected_h.pass_h_shapecorr-D3b.corrected_h.delta_h_true;temp(~isfinite(temp))=0;
D3b.ref_surf.h_corr_mean(sumW>0)=sum(W(sumW>0,:).*(temp(sumW>0,:)), 2);


% calculate the point-by-point errors
good_inc=D3b.pass_quality.summary==0 &  D3b.corrected_h.pass_included_in_fit;
good=D3b.pass_quality.summary==0;
%abs(D3b.pass_stats.y_RGT-D3b.corrected_h.ref_pt_y_RGT)<100;
dz_true=D3b.non_product.delta_z_firn+D3b.non_product.delta_h;
[sigma_dz, sigma_dz_inc, Ndz, Ndz_inc, sigma_dz_s, zbar]=deal(NaN(size(D3b.corrected_h.pass_h_shapecorr,1), 1));
for k=1:size(D3b.corrected_h.pass_h_shapecorr,1)
    igood=find(good(k,:));
    EDZ=(D3b.corrected_h.pass_h_shapecorr(k,good(k,:)))-dz_true(k, good(k,:));
    zbar(k)=interp1(igood, EDZ, 6);
    sigma_dz(k)=std(EDZ);
    ErrorEstSqr=D3b.corrected_h.pass_h_shapecorr_sigma(k,good(k,:)).^2+D3b.corrected_h.pass_h_shapecorr_sigma_systematic(k,good(k,:)).^2;
    sigma_dz_s(k)=sqrt(sum((EDZ-mean(EDZ)).^2./ErrorEstSqr)/(length(EDZ)-1));
    Ndz(k)=length(EDZ);
    EDZ=(D3b.corrected_h.pass_h_shapecorr(k,good_inc(k,:)))-dz_true(k, good_inc(k,:));
    sigma_dz_inc(k)=std(EDZ);
    Ndz_inc(k)=length(EDZ);
end

axes(hax1(2)); cla; hold on;
for k=1:12
    plot(D3b.reference_point.x_ATC(good(:,k)), D3b.corrected_h.pass_h_shapecorr(good(:,k),k)-zbar(good(:,k)),'.','color', colors(k,:)); 
end
set(gca,'ylim', [-14 14]);



axes(hax2(1)); cla; hold on;
[nn, bb]=histcounts(sigma_dz(isfinite(sigma_dz) & Ndz > 4), 0:.02:1);
hb(1)=bar(bb(1:end-1)+diff(bb)/2, nn,'k');
med_a=percentile(sigma_dz(isfinite(sigma_dz)), 0.68);
plot([1 1]*med_a, get(gca,'ylim')); 
ht=text(med_a, mean(get(gca,'ylim')), num2str(med_a));
 
[nn, bb]=histcounts(sigma_dz(isfinite(sigma_dz_inc) & Ndz_inc > 4), 0:.02:1);
hb(2)=bar(bb(1:end-1)+diff(bb)/2, nn,'r');
xlabel('error, m'); ylabel('count');

 
axes(hax2(2)); cla; 
histogram(sigma_dz_s(isfinite(sigma_dz_s)), 0:.05:3); set(gca,'xlim', [0 3]); hold on;
med_s=percentile(sigma_dz_s(isfinite(sigma_dz_s)), 0.68);
plot([1 1]*med_s, get(gca,'ylim')); 
xlabel('error/\sigma'); ylabel('count');
ht=text(med_s, mean(get(gca,'ylim')), num2str(med_s));

%plot(-6:.02:6, 0.02*sum(isfinite(temp(:)))*gaussian(-6:.02:6, 0, 1),'r','linewidth', 2)
 
linkaxes(hax1,'x');
set(findobj(hax1,'type','line'),'buttondownfcn', 'ATL06_11_bdf(D3)');
set( hax1(1),'buttondownfcn', 'ATL06_11_BDF(D3)');
