% Script - star_azel.m
% Computes look angles from input lat/lon to a star using data from
% the nautical almanac for given data and time.
%
% Inputs: 
%   lat -           Observer latitude (eg. 39.0 N)
%   lon -           Observer longitude (eg. 110.0 W)
%   SHA -           Sidereal hour angle (degrees west of vernal equinox)
%   dec -           declination (degrees north/south of celestial equator)
%   GHA_Aries -     Greenwich Mean Hour Angle for Aries
%
% Ouputs:
%   az -        azimuth to star from observer position
%   el -        elevation to star from observer position
%
% Author: drohm
%-------------------------------------------------------------------------
 
% Example: My location in Lone Tree, Colorado at 04:05:10(UTC) October
% 14th, 2020.
lat = 32.00;            % earth observer latitude
lon = -110.00;           % earth observer longitude  
yyyy = 2020;
mm = 12;
dd = 29;
HH = 19;
MM = 50;
SS = 00;

%-Example: Star is Arcturus, Convert data to decimal degrees 
SHA = 145 + 51.1/60;
dec = 19 + 04.4/60;
GHA_Aries = calc_GHA_Aries(yyyy,mm,dd,HH,MM,SS);

%-Compute Azimuth and Elevation angles from location to celestial object
GHA = GHA_Aries + SHA;
LHA = GHA + lon;   
%-Add or subtract multiples of 360 to get in range of 0 - 360
while LHA > 360
    LHA = LHA-360;
end
while LHA < 0 
    LHA = LHA + 360;
end

[Hc,Z] = calcHCZ(dec,lat,LHA);
Az = Z
El = Hc


