function y=percentile(x, p);

if length(x) < 2
    y=NaN; 
    return
end

s=sort(x); 
P0=(1:length(x))/length(x);
y=interp1(P0, s, p);