function out=platte_carre_WGS84(in, latlon0, xform)

WGS84.semimajor = 6378137;
WGS84.semiminor = 6356752.314245;
WGS84.flattening = 1/298.257223563;
WGS84.eccentricity = sqrt(0.00669437999013);
WGS84.ellipsoid = [WGS84.semimajor WGS84.eccentricity];

if ~exist('xform','var'); 
    xform='forward'; % forward, latlon->ne;
end

if ~ismember(xform,{'distance','track_azimuth','offset_by_az_dist'});
    lat0=latlon0(:,1);
    lon0=latlon0(:,2);
    scale_lat=(1-WGS84.eccentricity^2).*pi.*WGS84.semimajor./(180*((1-(WGS84.eccentricity*sind(lat0)).^2).^1.5));
    scale_lon=cosd(lat0).*pi.*WGS84.semimajor./(180*((1-(WGS84.eccentricity*sind(lat0)).^2).^0.5));
end

switch xform
    case 'forward'
        lat=in(:,1);
        lon=in(:,2);
        out(:,2)=(lat-lat0).*scale_lat;
        out(:,1)=(lon-lon0).*scale_lon;
    case 'reverse'
        x=in(:,1);
        y=in(:,2);
        out(:,1)=y./scale_lat+lat0;
        out(:,2)=x./scale_lon+lon0;
    case 'track_azimuth' % azimuth for latlon
        dxy1=platte_carre_WGS84(in(2:end,:), in(1:end-1,:),'forward');
        dxy2=-platte_carre_WGS84(in(1:end-1,:), in(2:end,:),'forward');
        dxy=(dxy1+dxy2)/2;  
        out=atan2d(dxy(:,1), dxy(:,2));
        out(end+1)=out(end);
    case 'distance'  % distance 
        dxy1=platte_carre_WGS84(in(2:end,:), in(1:end-1,:),'forward');
        dxy2=-platte_carre_WGS84(in(1:end-1,:), in(2:end,:),'forward');
        dxy=(dxy1+dxy2)/2;
        out=cumsum([0; sqrt(dxy(:,1).^2+dxy(:,2).^2)]);
    case 'offset_by_az_dist'
        az=latlon0(:,1);
        dist=latlon0(:,2);
        delta=[sind(az).*dist cosd(az).*dist];
        out=platte_carre_WGS84(delta, in,'reverse' );
        
end
        
