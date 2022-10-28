
%--------------------------------------------------
function ATL11_xtrack_plot(D3, hc_rep, Ymed, h0, m_surf, sigma_m_surf, fit_reps, params)

 axes(params.axes(1,2));

yvals=(min(D3.y_RGT):5:max(D3.y_RGT))';
Gtemp=[poly_fit_mat((yvals-Ymed)/params.y_scale, length(m_surf)-1)];
plot(yvals, h0+Gtemp*m_surf(2:end),'b');hold on;
surf_error=sqrt(sum((Gtemp*diag(sigma_m_surf(2:end))).^2,2));


plot(yvals, h0+surf_error,'b--')
plot(yvals, h0-surf_error,'b--')

plot(D3.y_RGT, D3.h_LI,'rx')
 
included=ismember(D3.rep, fit_reps);
plot(D3.y_RGT(included), D3.h_LI(included),'ko');
ylabel('height, m');

axes(params.axes(2,2));cla
plot(yvals,  Gtemp(:, 2:end)*m_surf(3:end),'b'); hold on;
plot(yvals,  Gtemp(:, 2:end)*m_surf(3:end)+surf_error,'b--')
plot(yvals,  Gtemp(:, 2:end)*m_surf(3:end)-surf_error,'b--')
dumbell_plot(D3.y_RGT, D3.h_LI-h0-(D3.y_RGT-Ymed)/params.y_scale*m_surf(2), D3.rep, D3.beam, 'linestyle','--');

dh_corr=hc_rep(D3.rep)';
dumbell_plot(D3.y_RGT, h0+D3.h_LI-h0-(D3.y_RGT-Ymed)/params.y_scale*m_surf(2)-dh_corr, D3.rep, D3.beam, 'linestyle',':');

dumbell_plot(D3.y_RGT(included), D3.h_LI(included)-h0-(D3.y_RGT(included)-Ymed)/params.y_scale*m_surf(2), D3.rep(included), D3.beam(included),'linestyle','--');
dumbell_plot(D3.y_RGT(included), h0+ D3.h_LI(included)-h0-(D3.y_RGT(included)-Ymed)/params.y_scale*m_surf(2)-dh_corr(included), D3.rep(included), D3.beam(included),'linestyle','--');
ylabel({'height' '(linear correction)'}')
