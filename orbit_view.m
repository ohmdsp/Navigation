function [sqmview, sqmlview, coverage, N] = orbit_view(height)

% Compute approximate satallite visibility in square miles base on orbit height.
%
% Inputs: 
%   height       -   Height of orbit (km)
%   
% Ouputs:
%   sqmview     -   Approximate visability in meters^2
%   sqmlview    -   Approximate visability in square miles
%   coverage    -   Approximate % coverage of earth by single satellite
%   N           -   Number of satellites needed to cover entire earth 
%
% Author: drohm
%-------------------------------------------------------------------------

%-Constants
Re = 6378.137;          % equitorial radius of the earth (km) - WGS84

%-Compute view area
theta = acos(Re/(Re + height));
s = Re*theta;
smiles = s / 1.60934;       % conversion from meters to miles 
sqmview = s^2;              % footprint in square km
sqmlview = smiles^2;        % footprint in square miles

%-Compute % coverage of earth surface by one satellite
R_area = 4*pi*Re^2;      % square km
coverage = sqmview/R_area * 100;     % percentage of earthcoverage 

%-Compute number of satellites needed to cover entire earth surface
N = ceil(100/coverage);
