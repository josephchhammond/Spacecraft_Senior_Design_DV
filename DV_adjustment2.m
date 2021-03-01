function [dtburn_dt] = DV_adjustment2(DV2_factor,p,dtburn_dt_max)
% This function adjusts DV based on non-instantaneous behavior
% 
% Inputs:
% DV2_factor - Capacity of system divided by needs of system
% p - Coefficients for adjustment equation for noninstantaneous behavior
% dtburn_dt - Max fraction of transfer spent in burn

% Output:
% dtburn_dt - Minimum fraction of transfer spent in burn

%DV_adj = DV * (polyval(p,dtburn_dt)).^.5; %original
a = p(1);
b = p(2);
c = p(3) - DV2_factor^2;

[x1,x2] = quad(a,b,c);

check1 = zeros(size(DV2_factor));
check2 = dtburn_dt_max*ones(size(DV2_factor));

% Find positive solution
check_x1 = x1 > check1;
check_x2 = x2 > check1;
dtburn_dt = x1.*check_x1 + x2.*check_x2;

% Set any values above dtburn_dt_max to dtburn_dt_max
excess_check = dtburn_dt > check2;
dtburn_dt = dtburn_dt - excess_check.*dtburn_dt + dtburn_dt_max.*excess_check;
end
