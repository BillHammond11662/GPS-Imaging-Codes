function [delta,az]=baz(lat1,lon1,lat2,lon2)
%  function [delta,az]=baz(lat1,lon1,lat2,lon2)
%
%  Calculates the delta and azimuth between two points on Earth
%
%  all inputs are in degrees
%
%  performs similar function to "distance.m" in Matlab, with much
%  faster run time.

dtor = pi/180;
lat1 = dtor*lat1;
lat2 = dtor*lat2;
lon1 = -dtor*lon1;
lon2 = -dtor*lon2;

cosdla = cos(lon2-lon1);

ss = sin(lat1).*sin(lat2);
cc = cos(lat1).*cos(lat2);
delta = acos(ss + cc.*cosdla);

ss = sin(lat2)-sin(lat1).*cos(delta);
cc = sin(delta).*cos(lat1);

az = acos(ss./cc);

delta = delta./dtor;
delta = real(delta);

az = az./dtor;

if (sin(lon2-lon1) > 0)
    az = 360 - az;
end

if (delta>180 | delta<0)
    error ('problems arose');
end
