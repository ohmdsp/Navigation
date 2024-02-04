function [dist_km, dist_mile] = great_circle(lat1, lon1, lat2, lon2)

% Computes the great circle distance between two lat/lon positions on
% surface of the earth.
%
% Inputs: 
%   lat1 -          position 1 latitude 
%   lon1 -          position 1 longitude
%   lat2 -          position 2 latitude 
%   lon2 -          position 2 longitude
%
% Ouputs:
%   dist_km -       great circle distance in kilometers
%   dist_mile -     great circle distance in miles 
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

%-Compute dot product and calulate distance miles
%theta_c = acos(dot(A,B)./(norm(A)*norm(B)));
theta_c = acos(dot(A,B)./(R^2));

dist_km = R*theta_c;
dist_mile = R*theta_c * .621371;
