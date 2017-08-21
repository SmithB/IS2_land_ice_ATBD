function count=my_histc(x, bins)
x=sort(x);
bins=bins(:);
if all(isnan(bins)); count=NaN; return; end

% case 1: bins are of equal spacing
% NEW: October 2015: Better check on bin width
 if abs(max(diff(bins))-min(diff(bins)))/mean(diff(bins)) < 1e-9
    ind=round((x-bins(1))/(bins(2)-bins(1)))+1;
    [bin_num, ia]=unique(ind,'first');
    [~, ib]=unique(ind,'last');
    count=zeros(size(bins));
    % NEW: October 2015: handle single-event bins
    good=isfinite(bin_num) & bin_num > 0 & bin_num <= length(count);
    count(bin_num(good))=ib(good)-ia(good)+1;
else
    count=zeros(size(bins));
    bin_edges=[bins(1)-(bins(2)-bins(1))/2; 
        (bins(2:end)+bins(1:end-1))/2; 
        bins(end)+(bins(end)-bins(end-1))/2];
    for k=1:length(bin_edges)-1;
        count(k)=sum(x>=bin_edges(k) & x <bin_edges(k+1));
    end
end
    
