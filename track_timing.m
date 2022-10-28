delta_t_orb=(datenum('14-Oct-2018 23:12:26')-datenum('14-Oct-2018 01:12:24'))/(264-250);


t0=datenum('16-Oct-2018 16:22:46');
orb0=280;
orb_num=orb0:orb0+500;



for k=1:length(orb_num)
    fprintf(1, '%s %d\n', datestr(t0+delta_t_orb*(orb_num(k)-orb0)), orb_num(k)); 
end


