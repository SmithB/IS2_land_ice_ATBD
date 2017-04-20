function [x0, search_hist]=linesearch_q5(f, XR, params)

if ~exist('params','var')
    params.xtol=diff(range(XR))/1e4;
end

if length(XR)==2
    xx=linspace(XR(1), XR(2), 10);
else
    xx=XR;
end

search_hist.x=[]; search_hist.z=[];

count=0;
while diff(range(xx)) > params.xtol  || count==0
    count=count+1;
    yy=f(xx);    
    best=find(yy==min(yy), 1, 'first');
    if best==1; best=2; end; 
    if best==length(xx); best=length(xx)-1; end
    search_hist.x=[search_hist.x; xx(:)];
    search_hist.z=[search_hist.z; yy(:)];
    xx_last=xx;
    xx=linspace(xx(best-1), xx(best+1), 10);
end
x0=xx_last(best);
