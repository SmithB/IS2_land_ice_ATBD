function z=dilate_mask(z, width);

%if(width/2==floor(width/2));
    [gx, gy]=meshgrid(-width/2+1/2:width/2-1/2);
%else
%    [gx, gy]=meshgrid(-width/2:width/2);
%end
kernel=double(abs(gx+i*gy)<=width/2);


z=conv2(double(z), kernel, 'same')>0;

