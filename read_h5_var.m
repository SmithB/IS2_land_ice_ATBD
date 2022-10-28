%--------------------------------
function x=read_h5_var(thefile, thevar, start, count);


if exist('start','var');
    x=double(h5read(thefile, thevar, start, count));
else
    x=double(h5read(thefile, thevar));
end
try
    temp=double(h5readatt(thefile,thevar,'_FillValue'));
    x(x==temp)=NaN;
end

