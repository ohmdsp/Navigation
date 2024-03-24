function [range_km, range_mile, bearing_deg] = range_and_bearing(lat1, lon1, lat2, lon2)

% Computes the range and bearing between two lat/lon positions on
% surface of the earth using (1) spherical geometry solution. 
% Note: within 100nmi (1nmi = 1580meters), can assume flat earth with 
% error less than 0.3nmi.
%
% Inputs: 
%   lat1 -          position 1 latitude 
%   lon1 -          position 1 longitude
%   lat2 -          position 2 latitude 
%   lon2 -          position 2 longitude
%
% Ouputs:
%   dist_km -       range in kilometers
%   bearing -       initial bearing along great circle measured clockwise
%                   in degrees from north
%
% Author: drohm
%-------------------------------------------------------------------------

%-Constants
R = 6378.137;                % equitorial radius of the earth (km) - WGS84

%-Convert positions to ECEF coordinates
az = R*sind(lat1);
ax = R*cosd(lat1)*cosd(lon1);
ay = R*cosd(lat1)*sind(lon1);

bz = R*sind(lat2);
bx = R*cosd(lat2)*cosd(lon2);
by = R*cosd(lat2)*sind(lon2);

%-Form vectors from ECEF components
A = [ax ay az]';
B = [bx by bz]';

%-Compute dot product and calulate range
%theta_c = acos(dot(A,B)./(norm(A)*norm(B)));
%theta_c = acos(dot(A,B)./(R^2));
theta_c = acos(sind(lat1)*sind(lat2)+cosd(lat1)*cosd(lat2)*cosd(lon1-lon2));

range_km = R*theta_c;
range_mile = R*theta_c * .621371;

%-Compute bearing
y = sind(lon2-lon1)*cosd(lat2);
x = cosd(lat1)*sind(lat2) - sind(lat1)*cosd(lat2)*cosd(lon2-lon1);
bearing_tmp = atan2d(y,x);
bearing_deg = mod(bearing_tmp,360);
