function [Az,El,rho] = geostat_azel(theta_e, phi_e, phi_s, hemisphere)

% Computes look angles from input lat/lon to geostationay satellite.
%
% Inputs: 
%   theta_e -       Observer latitude (eg. 39.0 N)
%   phi_e -         Observer longitude (eg. 110.0 W)
%   phi_s -         satellite longitude (eg. 119 for DirectTV 7S) 
%   hemisphere -    1 = northern hemisphere, 2 = southern hemisphere
%
% Ouputs:
%   az -        Satellite azimuth from observer position
%   el -        Satellite elevation from observer position
%   rho -       range to satellite from observer position
%
% Author: drohm
%-------------------------------------------------------------------------

%-Constants
R = 6378.137;                % equitorial radius of the earth (km) - WGS84
rs = 42200;                  % mean radius of geostationary orbit (km)

%-Compute difference in satellite and observer longitude
phi_se = phi_s - phi_e;    

%-Compute satellite azimuth and elevation from specified earth station using
%-Richharia's equations (Satellite Communication Systems, 2nd ed., M. Richharia) 
delta = 0;                  % declination (or inclination)
eta_s = asind( sind(delta)*sind(theta_e) + cosd(delta)*cosd(theta_e)*cosd(phi_se) );
El = atand( (sind(eta_s) - R/rs)/cosd(eta_s) );
Az = atand( sind(phi_se)/(cosd(theta_e)*tand(delta) - sind(theta_e)*cosd(phi_se)) );

%-Make corrections for hemisphere and earth azimuth quadrant
if hemisphere == 1          % northern hemisphere
    if phi_se < 0
        Az = 180 + Az;
    elseif phi_se > 0
        Az = 180 - abs(Az);
    else
        Az = Az;
    end   
elseif hemisphere == 2      % southern hemisphere
    if phi_se < 0
        Az = 360 - Az;
    elseif phi_se > 0
        Az = 180 - Az;
    else
        Az = Az;
    end
end
   
%-Compute range from observer to satellite
rho = sqrt( rs^2 - R^2*cosd(El)^2 ) - R*sind(El);