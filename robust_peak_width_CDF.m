function [sigma_hat, med]=robust_peak_width_CDF(x, N_BG, XR)

% C_tot=C_sig+C_noise
%      =C_sig(t)+(t-t0)BGR
% want t1 such that C_sig(t1)=f
% Ctot(t1)=f+(t1-t0)BGR
% look for first point where Ctot > 0.16 Nsig + (t1-t0)*BGR 
% and last point where Ctot < 0.84 Nsig + (t1-t0)*BGR
% the difference gives the peak width

% if ~exist('XR','var')
%     XR=[min(x) max(x)];
% end
BGR=N_BG/diff(XR);

N_sig=length(x)-N_BG;

xs=sort(x); xs=xs(:);
C=0.5:length(x)-0.5; C=C(:);

%new October 2015: switch to P vals of .25 and .75, include correction
%based on the erfinv 
 P_vals=[.25 .75];
% scale_factor=diff(erfinv(P_vals));
% N.B.  This is always the same, so compute by hand:
scale_factor=.5881;

i0=find(C<P_vals(1)*N_sig + (xs-XR(1))*BGR, 1, 'last'); 
if isempty(i0); i0=1;end
i1=find(C>P_vals(2)*N_sig + (xs-XR(1))*BGR, 1, 'first');
if isempty(i1); i1=length(x);end
sigma_hat=diff(xs([i0 i1]))/2/scale_factor;

if nargout==2;
    i0=find(C<0.5*N_sig + (xs-XR(1))*BGR, 1, 'last');
    if isempty(i0); i0=1;end
    i1=find(C>0.5*N_sig + (xs-XR(1))*BGR, 1, 'first');
    if isempty(i1); i1=length(x);end
    
    med=interp1(C([i0 i1])-(xs([i0 i1])-XR(1))*BGR, xs([i0 i1]), 0.5*N_sig);
end

