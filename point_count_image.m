function [count, xg, yg]=point_count_image(x, y, x0, y0)


if length(x0)==1
    xr=round_to(range(x), x0);
    x0=[xr(1)-x0 x0 xr(2)+x0];
end
if length(y0)==1
   yr=round_to(range(y), y0);
   y0=[yr(1)-y0 y0 yr(2)+y0];
end


ind_x=ceil((x-x0(1))/x0(2));
ind_y=ceil((y-y0(1))/y0(2));

xg=x0(2)+(x0(1):x0(2):x0(3)); 
yg=y0(2)+(y0(1):y0(2):y0(3)); 

nx=length(xg);
ny=length(yg);

good=ind_x>0 & ind_x <=nx & ind_y>0 & ind_y < ny;
if any(~good)
    ind_x=ind_x(good); 
    ind_y=ind_y(good);
end

count=sparse([], [], [], ny, nx);
ind_bins=0:1e6:length(ind_y);
if ind_bins(end) < length(ind_y); 
    ind_bins(end+1)=length(ind_y);
end
for k=1:length(ind_bins)-1;
    these_bins=ind_bins(k)+1:ind_bins(k+1);
    count=count+sparse(ind_y(these_bins), ind_x(these_bins), ones(length(these_bins),1), ny, nx);
end

