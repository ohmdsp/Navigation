function [Hc,Z] = calcHCZ(dec,lat, LHA)

% Calculates Hc (elevation) and Z (azimuth) for celestial navigation
% Inputs assumed to be in decimal degrees
%
% Inputs: 
%   dec -       declination
%   lat -       Observer latitude (degrees)
%   LHA -       Local Hour Angle
%
% Ouputs:
%   Hc -        Star elevation from observer position
%   Z -         Star azimuth relative to observer position
%
% Author: drohm
%-------------------------------------------------------------------------

S = sind(dec);
C = cosd(dec)*cosd(LHA);
Hc = asind( S*sind(lat)+C*cosd(lat) );

X = (S*cosd(lat) - C*sind(lat))/cosd(Hc);
if X>1
    X = 1;
elseif X<-1
    X = -1;
else
    X = X;
end

A = acosd(X);

if LHA>180
    Z = A;
else
    Z = 360-A;
end


end

