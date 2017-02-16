function varargout=polar_stereo_xform(in1, in2, xform, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% turns lat,lon pairs into polar stereographic pairs
% uses E=.082271850
% R=6378.2650 km
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


%constants:
PI  =	4*atan(1);
DR  = PI/180;
CDR = 180/PI;
RE=6378137.0;  % for WGS84
E   = 0.08181881069348491;% for WGS84; UPDATED 11/2012
E2 = 0.00669431778329633;
SLAT=abs(params.standard_latitude);
SN=params.hemisphere;
lon_offset=params.lon_offset;

t   = tan(PI/4.0 - SLAT/(2.0*CDR)) ./ ( ((1.0 - E*sin(DR*SLAT))/ (1.0 + E*sin(DR*SLAT)))^(E/2.0) );
cm  =	cos(DR*SLAT)/sqrt(1.0 - E2*(sin(DR*SLAT)*sin(DR*SLAT)));
 
if strcmp(xform,'forward')
    lat=in1;
    lon=in2;
    t1 = tan ( PI / 4.0 - SN*lat / ( 2.0 * CDR) )./( ( ( 1.0 - E * sin( DR * SN*lat ) ) ./ ( 1.0 + E * sin( SN*DR * lat ) ) ) .^ ( E / 2.0 ));
    rho = RE * cm * t1 / t;
    lon_offset=params.lon_offset*DR;
    x = real(rho.*sin(DR*lon+lon_offset));
    y = real(rho.*cos(DR*lon+lon_offset));

    varargout={x, y};

elseif strcmp(xform,'reverse');
    x=in1; 
    y=in2;
    
    t 	= 	tan(PI/4.0 - SLAT/(2.0*CDR)) / ( ((1.0 - E*sin(DR*SLAT))/(1.0 + E*sin(DR*SLAT)))^(E/2.0) );
    cm 	= 	cos(DR*SLAT)/sqrt(1.0 - E2*(sin(DR*SLAT)*sin(DR*SLAT)));
    SN=params.hemisphere;
    rho = sqrt( x .* x + y .* y );
    chi = PI/2.0 - 2.0 * atan2( t * rho, RE * cm );
    lat = chi +...
        ((E2/2.0) + (5.0*E2*E2/24.0) + (E2*E2*E2/12.0))*sin(2.0*chi) + ...
        ((7.0*E2*E2/48.0) + (29.0*E2*E2*E2/240.0))*sin(4.0*chi) + ...
        (7.0*E2*E2*E2/120.0)*sin(6.0*chi);
    lat = lat*(SN * CDR);
    lon = params.lon_offset + CDR*atan2(x, y);
    lon((lon < -180))= lon(lon<-180)+360.0;
    lon(lon > 180.0 )=lon(lon>180.0) - 360.0;
    varargout={lat, lon};
else
    error(sprintf('xform must be ''forward'' or ''reverse.'' ''%s'' not understood.', xform));
end
