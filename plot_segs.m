function h=plot_segs(x, h, slope, W, linespec);


plot([x(:)-W/2 x(:)+W/2]', [h(:)-slope(:)*W/2 h(:)+slope(:)*W/2]', linespec);