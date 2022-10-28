function hh=plot_segs(x, h, slope, W, linespec)

xx=[x(:)-W/2, x(:)+W/2, nan(size(x(:)))].';
yy=[h(:)-slope(:)*W/2, h(:)+slope(:)*W/2, nan(size(h(:)))].';
hh=plot(xx(:), yy(:), linespec);
%h=plot([x(:)-W/2 x(:)+W/2]', [h(:)-slope(:)*W/2 h(:)+slope(:)*W/2]', linespec);