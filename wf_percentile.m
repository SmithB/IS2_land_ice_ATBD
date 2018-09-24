function X=wf_percentile( bins, counts, P)

if isempty(bins); X=NaN; return; end;

bin_width=[diff(bins(:)); bins(end)-bins(end-1)];
edges=[bins(1)-bin_width(1)/2; bins(:)+bin_width/2];
C=[0; cumsum(counts(:))]; C=C/C(end);
for k=1:length(P)
    % check if any C are equal to P(k)
    eq_els=abs(C-P(k)) < 100*eps;
    if any(eq_els)
        X(k)=mean(edges(eq_els));
        continue
    end
    i_minus=find(C<P(k)-100*eps, 1, 'last');
    i_plus=find(C>=P(k)+100*eps, 1, 'first');
    if C(i_plus)==C(i_minus)
        X(k)=0.5*(edges(i_plus)+edges(i_minus));
    else
        X(k)=((P(k)-C(i_minus))*edges(i_plus)+(C(i_plus)-P(k))*edges(i_minus))/(C(i_plus)-C(i_minus));
    end
end



 