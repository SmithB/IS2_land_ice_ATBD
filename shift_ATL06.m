function D3a=shift_ATL06(D3a,TrackData, shift_mag)

az0=interp1(TrackData.x_RGT, TrackData.azimuth, D3a.x_RGT(:,1),'pchip');

for k=1:2;
    ll=platte_carre_WGS84(  [D3a.lat_ctr(:,k), D3a.lon_ctr(:,k)], [az0-90, shift_mag*ps_scale_for_ll(D3a.lat_ctr(:,k))],'offset_by_az_dist');
    D3a.lat_ctr(:,k)=ll(:,1); 
    D3a.lon_ctr(:,k)=ll(:,2); 
end
xx=ll2ps(D3a.lat_ctr, D3a.lon_ctr);

D3a.x_PS_ctr=real(xx); 
D3a.y_PS_ctr=imag(xx);
D3a.y_RGT=D3a.y_RGT+shift_mag;