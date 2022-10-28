function varargout=ll2ps(lat, lon)

% Antarctic polar-stereographic forward transform, with standard latitude of -71,
% central longitude of 0.  This is a wapper function that calls
% 'polar_stereo_xform.m' with the appropriate parameters


[x,y]=polar_stereo_xform(lat, lon,'forward', struct('standard_latitude', -71,'hemisphere', -1, 'lon_offset', 0));
if nargout==1
    varargout={x+1i*y};
elseif nargout==2
    varargout={x y};
else
    error('ll2ps:too many output arguments');
end