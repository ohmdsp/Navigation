function GHA_Aries = calc_GHA_Aries(yyyy,mm,dd,HH,MM,SS)

% Calculates GHA_Aries

% Inputs:
%   yyyy -          Year of observation
%   mm -            Month of observation (1-12)
%   dd -            Day of observation (numerical value)
%   HH -            Hour of observation (0-24)
%   MM -            Minute of observation (0-60)
%   SS -            Second of observation (0-60)

% To calculate the GHA for Aries, first calculate the Julian Day, 
% which is the number of elapsed days (and fractions of a day) since 
% beginning of the Julian Period (12 am, 1 January, 4713 BC) which is a 
% common reference epoch in celestial navigation. 
% The Julian Day for yyyy-mm-dd HH:MM:SS is calculated using the formula:
JD = (367 * yyyy - floor(7*(yyyy + floor((mm+9)/12))/4) - ...
    floor(3*(floor((yyyy+(mm-9)/7)/100)+1)/4) + floor(275*mm/9) + ...
    dd + 1721028.5) + (HH+MM/60+SS/3600) / 24;

% Next, we calculate a factor (T) referred to as Julian Centuries which is
% a unit of elapsed time since 1 January, 2000, 12 am. The formula is:
T = (JD - 2451545 ) / 36525;

% Now we can calculate the GHA for Aries at the time of our observation
GHA_Aries = mod( (280.46061837 + 360.98564736629 * (JD-2451545) + ...
    0.000387933*T*T - T*T*T/38710000), 360);

% The formula calculates the angular “movement” of Aries along the 
% celestial sphere since the established epoch (J2000, i.e. 1 January 2000, 12 am)
% and the modulus of 360 brings the value to an angle in degrees.

