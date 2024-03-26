function [height, radius] = orbit_height(T)

% Compute satallite orbit height and radius using pendulum method. 
%
% Inputs: 
%   T       -   Period of orbit (seconds)
%   
% Ouputs:
%   height  -   Satellite height above the surface of earth (meters)
%   radius   -   Satellite orbit radius (meters)
%
% Author: drohm
%-------------------------------------------------------------------------

%-Constants
Re = 6378.137;          % equitorial radius of the earth (km) - WGS84
G = 6.674e-11;          % Gravitational constant (m^3/(kg*s^2))
Me = 5.972e24;          % mass of the earth (kg)

%-Compute orbit radius (km) and height (km) using pendulum method
radius = ( G*Me*T^2*(1/(2*pi))^2 )^(1/3) / 1000;
height = radius - Re;
