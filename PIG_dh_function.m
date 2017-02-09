function dh=PIG_dh_function(x, y, days)

lambda=30000;
amp=10;

xy_fun=sin(2*pi*x/lambda).*sin(2*pi*y/lambda); 

t_fun=zeros(size(days));
years=days/365;
these=(years > 1 & years <=2);
t_fun(these)=years(these)-1;
t_fun(years>2)=1;


dh=amp.*xy_fun.*t_fun;
 