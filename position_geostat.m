% position_geostat.m
% Find position from intercept and azimuth for two or more
% geostationary satellites. See 2020 Nautical Almanac pg. 282.
%
% Author: David Ohm
% -----------------------------------------------------
clear all;close all;clc

%-Constants
R = 6378.137;                % equitorial radius of the earth (km) - WGS84
rs = 42200;                  % mean radius of geostationary orbit (km)

%-True position
lat = 39.1427;
lon = -76.8606;

%-Azimuth and Elevation Observations
[Zo1,Ho1,rho1] = geostat_azel(lat, lon, -115.0, 1);
[Zo2,Ho2,rho2] = geostat_azel(lat, lon, -27.5, 1);
[Zo3,Ho3,rho3] = geostat_azel(lat, lon, -50.0, 1);

%-Starting coordinates (Assumed Position)
Bf = 41;                % observar latitude  
Lf = -75;               % observer longitude

%-Start iteration
d = 100;
while d > 20

    %-Intercept and azimuth calculations (based on assumed position)
    [Z1,Hc1,rho1] = geostat_azel(Bf, Lf, -115, 1);
    [Z2,Hc2,rho2] = geostat_azel(Bf, Lf, -27.5, 1);
    [Z3,Hc3,rho3] = geostat_azel(Bf, Lf, -50.0, 1);
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

    Li = Lf + (A*E - B*D)/(G*cosd(Bf));
    Bi = Bf + (C*D - B*E)/G;

    %-Calculate distand between initial estimated position and new estimate
    d = 60*sqrt( ((Li-Lf)^2*cosd(Lf)^2 + (Bi - Bf)^2 ) );
    
    Lf = Li;
    Bf = Bi;
    
end

disp(['Estimated position is latitude = ' num2str(Bf) ', longitude = ' num2str(Lf) ])


  

