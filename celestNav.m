% Script clestNav.m
% Algorithm for computing position estimated by celestial navigation
% (checked using 2018 and 2020 nautical almanac example)
%
% Author: drohm
%------------------------------------------------------------------------
clear all;close all;clc

%-Date and time (UT) of last position fix and lat/lon (degrees)
% 2020 July 6, 21:00:00 UT
lat = 32;
lon = -15;

%-Stars with corrected (observed) altitudes and observation times
% (1)Regulus -> 20:39:23 -> 24.9810
% (2)Antares -> 20:45:47 -> 26.8969
% (3)Kochab -> 21:10:34 -> 47.4869
Ho1 = 24.9810;
Ho2 = 26.8969;
Ho3 = 47.4869;

%-Ship speed (knots) and heading (degrees)
v_s = 20;
v_h = 325;

%-GHA for celestial objects interpolated to time of observation
%   Look up the GHA-Aires for the integer hour prior and hour after time of
%   observation. Interpolate to time of observation using
%   GHA = GHA_0 + x*(GHA_1 - GHA_0) where x = min/60+sec/60 of time past
%   integer hour prior to observation (x will be less than 1).
%
%   For example, 20h GHA_0 = 210deg19.0 from nautical almanac. This converts
%   210.3167 decimal degrees. Likewise 21h GHA_1 = 225.3583. If
%   obeservation of Regulus was made at time is at 20:39:23, so 
%   x = 39/60+23/360 = 0.6564. Therefore, GHA = 220.1840
%
%   Then compute GHA1 (GHA for Regulus) by GHA1 = GHA+SHA. Look up SHA 
%   for Regulus in Nautical Almanac.
yyyy = 2020;
mm = 7;
dd = 6;
HH = 20;
MM = 39;
SS = 23;
GHA_Aries1 = calc_GHA_Aries(yyyy,mm,dd,HH,MM,SS);

yyyy = 2020;
mm = 7;
dd = 6;
HH = 20;
MM = 45;
SS = 47;
GHA_Aries2 = calc_GHA_Aries(yyyy,mm,dd,HH,MM,SS);

yyyy = 2020;
mm = 7;
dd = 6;
HH = 21;
MM = 10;
SS = 34;
GHA_Aries3 = calc_GHA_Aries(yyyy,mm,dd,HH,MM,SS);


SHA1 = 207+38.6/60;
SHA2 = 112+20.0/60;
SHA3 = 137+19.4/60;

GHA1 = mod(GHA_Aries1 + SHA1,360);
GHA2 = mod(GHA_Aries2 + SHA2,360);
GHA3 = mod(GHA_Aries3 + SHA3,360);

% GHA1 = 82.7715;
% GHA2 = 349.0660;
% GHA3 = 20.2687;

%-Declination for each star interpolated to time of observation.
%   Look up dec for star in Nautical Almanac for observation date and
%   convert to decial degrees.
dec1 = 11.8700;
dec2 = -26.4767;
dec3 = 74.0783;

%-Convert time to decimal seconds, the compute lat/lon of ship at time 
% of observations given speed and heading
t0 = 21+0/60+0/3600;
t1 = (20+39/60+23/3600)-t0;
t2 = (20+45/60+47/3600)-t0;
t3 = (21 + 10/60 + 34/3600)-t0;

tmp = 100;
while tmp>20      % iterate until error is less than 20 nautical miles

    lon1 = lon + t1*(v_s/60)*sind(v_h)/cosd(lat);
    lon2 = lon + t2*(v_s/60)*sind(v_h)/cosd(lat);
    lon3 = lon + t3*(v_s/60)*sind(v_h)/cosd(lat);
    lat1 = lat + t1*(v_s/60)*cosd(v_h);
    lat2 = lat + t2*(v_s/60)*cosd(v_h);
    lat3 = lat + t3*(v_s/60)*cosd(v_h);

    %-LHA for stars
    LHA1 = GHA1 + lon1;
    LHA2 = GHA2 + lon2;
    LHA3 = GHA3 + lon3;

    %-Call function to compute predicted elevation and azimuth
    [Hc1,Z1] = calcHCZ(dec1, lat1, LHA1);
    [Hc2,Z2] = calcHCZ(dec2, lat2, LHA2);
    [Hc3,Z3] = calcHCZ(dec3, lat3, LHA3);

    %-Position line intercept
    p1 = Ho1 - Hc1;
    p2 = Ho2 - Hc2;
    p3 = Ho3 - Hc3;

    %-Position calulation from intercept and azimuth
    A = cosd(Z1)^2 + cosd(Z2)^2 + cosd(Z3)^2;
    B = cosd(Z1)*sind(Z1) + cosd(Z2)*sind(Z2) + cosd(Z3)*sind(Z3);
    C = sind(Z1)^2 + sind(Z2)^2 + sind(Z3)^2;
    D = p1*cosd(Z1) + p2*cosd(Z2) + p3*cosd(Z3);
    E = p1*sind(Z1) + p2*sind(Z2) + p3*sind(Z3);
    G = A*C - B^2;

    Lt = lon + (A*E - B*D)/(G*cosd(lat));
    Bt = lat + (C*D - B*E)/G;

    %-Calculate distand between initial estimated position and new estimate
    d = 60*sqrt( ((Lt-lon)^2*cosd(lat)^2 + (Bt - lat)^2 ) );
    tmp = d;                 % repeat using least squares if d>20nm
    lon = Lt;
    lat = Bt;
    
end

disp(['Estimated position is lattitude = ' num2str(lat) ', longitude = ' num2str(lon) ])



