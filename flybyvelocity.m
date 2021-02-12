function [v] = flybyvelocity(R,p)
% Flyby velocity for each heliocentric range is determined
% 
% Inputs:
% R = heliocentric range, AU
% p = polynomial coefficients describing flyby velocity
% 
% Outputs:
% v = flyby velocity, km/s
% 
% Written by: Joey Hammond

v = polyval(p,1/R);
end

