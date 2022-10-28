function X=wf_percentile_bnds( bins, counts, P, bnds)

if isempty(bins); X=NaN; return; end;

% E01  E12   E23     E34      E45  E56 E67
% |    |     |       |        |    |   |
%   B1    B2     B3       B4    B5   B6 
% 0    B1   B1+B2                     sum(B)

counts=counts(:); 
bins=bins(:);

bin_width=[diff(bins(:)); bins(end)-bins(end-1)];
%edges=[bins(:)-bin_width(:)/2; bins(end)+bin_width(end)/2];
edges=[bins(:) - bin_width(:)/2 , bins(:) + bin_width(:)/2];

C=[0; cumsum(counts(:))]; 
C=C/C(end);

if exist('bnds','var')
    
    % get the percentiles for the bounds
    Cbnds=interp1([edges(1,1); edges(:,2)], C, bnds);  
     
    if bnds(1) <= edges(1) 
        Bfirst = 1;
        Cfirst=0; 
    else
        Bfirst=find(edges(:,2) > bnds(1), 1, 'first');
        Cfirst=Cbnds(1);
    end
    
     
    if bnds(2) >= edges(end)  
        Blast=length(bins); Clast=1;
    else
        Blast=find(edges(:,1) < bnds(2), 1,'last');
        Clast=Cbnds(2);
    end
    
    edges=[bnds(1) edges(Bfirst,2); edges(Bfirst+1:Blast-1,:); edges(Blast,1), bnds(2)];
    C_edges=[[Cfirst; C(Bfirst+1:Blast)], [C(Bfirst+1:Blast); Clast]];
    C_edges=(C_edges-C_edges(1,1))/(C_edges(end,2)-C_edges(1,1));
end

X=zeros(size(P));

for k=1:length(P)
    % check if any C are equal to P(k)
    eq_els=abs(C_edges-P(k)) < 100*eps;
    if any(eq_els)
        X(k)=mean(edges(eq_els));
        continue
    end
    i_minus=find(C_edges(:,1)<P(k)-100*eps, 1, 'last');
    i_plus=find(C_edges(:,2)>=P(k)+100*eps, 1, 'first');
    if C_edges(i_plus,2)==C_edges(i_minus,1)
        X(k)=0.5*(edges(i_plus)+edges(i_minus));
    else
        X(k)=((P(k)-C_edges(i_minus,1))*edges(i_plus,2)+(C_edges(i_plus,2)-P(k))*edges(i_minus,1))/(C_edges(i_plus,2)-C_edges(i_minus,1));
    end
end



 