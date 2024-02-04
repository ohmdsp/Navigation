% position_geostat_newtons.m
%
% Find position from azimuth and elevation observations for two or more
% geostationary satellites using Newtons Method.
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

phi_s1 = -115;                  % sat1 longitude
phi_s2 = -27.5;                 % sat2 longitude
elev_1 = 30.1;                  % sat1 elevation measured
elev_2 = 22.3;                  % sat2 elevation measured

%-Starting coordinates (Assumed Position)
theta_e = 40;             % observar assumed latitude  
phi_e = -76;              % observer assumed longitude

%-Measured eta (left side of Riccharia's equation B.45)
meas1 = tand(elev_1);
meas2 = tand(elev_2);


% Start iteration
for i = 1:100
    
    X = [theta_e, phi_e];
    
    phi_se1 = phi_s1 - phi_e;
    phi_se2 = phi_s2 - phi_e;
    
    %-Note: geo satellite declination = 0
    eta_s1 = asind( cosd(theta_e)*cosd(phi_se1) ); % B.46
    eta_s2 = asind( cosd(theta_e)*cosd(phi_se2) );
    
    %-Compute predicted elevation
    eta1 = atand( sind(eta_s1) - R/rs )./ (cosd(eta_s1) );
    eta2 = atand( sind(eta_s2) - R/rs )./ (cosd(eta_s2) );
    
    %-df/dtheta_e
%     K11 = (sind(dec1)*cosd(theta_e) - cosd(dec1)*sin(theta_e)*sin(phi_se1))./cosd(eta_s1);    
%     K12 = (sind(dec2)*cosd(theta_e) - cosd(dec2)*sin(theta_e)*sin(phi_se2))./cosd(eta_s2); 
%     K21 = (cosd(dec1)*cosd(theta_e)*sind(phi_se1))./cosd(eta_s1);
%     K22 = (cosd(dec2)*cosd(theta_e)*sind(phi_se2))./cosd(eta_s2);
    K11 = (sin(theta_e)*sin(phi_se1))./cosd(eta_s1);    
    K12 = (sin(theta_e)*sin(phi_se2))./cosd(eta_s2); 
    K21 = (cosd(theta_e)*sind(phi_se1))./cosd(eta_s1);
    K22 = (cosd(theta_e)*sind(phi_se2))./cosd(eta_s2);
  
    %-Compute Jacobian Matrix
    A11 = K11*(1+tand(eta1)*tand(eta_s1) );
    A12 = -1*K21*(1-tand(eta1)*tand(eta_s1) );
    A21 = K12*(1+tand(eta2)*tand(eta_s2) );
    A22 = -1*K22*(1-tand(eta2)*tand(eta_s2) );
    J = [A11 A12;A21 A22];
    
    % Compute objective functions
    F1 = tand(eta1) - meas1;
    F2 = tand(eta2) - meas2;

    % Compute Xnew
    F = -1*[F1 F2];
    Jinv = inv(J);
    DeltaX = Jinv*F';
    Xnew = X+DeltaX'
    
    theta_e = Xnew(1);
    phi_e = Xnew(2);
    
end 


