function [dtburn_dt,DV_factor_adj] = DV_adjustment2(DV_factor,p,dtburn_dt_max)
% This function finds the burn/transfer time fraction for a given
% non-instantaneous burn based on a DV_factor, polynomial, and limit.
% Additionally, if the limit of burn/transfer time fraction is reached and
% adjusted DV_factor based on this limit is found.
% 
% Inputs:
% DV2_factor - Capacity of system divided by needs of system
% p - Coefficients for adjustment equation for noninstantaneous behavior
% dtburn_dt - Max fraction of transfer spent in burn

% Output:
% dtburn_dt - Minimum fraction of transfer spent in burn

% Original equation follows this form:
% DV_factor = (polyval(p,dtburn_dt)).^.5;

% Which can be rearranged as:
% polyval(p,dtburn_dt) - DV_factor^2 = 0;
% or
% a*dtburn_dt^2 + b*dtburn_dt + c - DV_factor^2 = 0




% Find quadratic solutions of dtburn_dt for a given DV_factor
a = p(2);
b = p(3);
c = p(4) - DV_factor.^2;
[x1,x2] = quad(a,b,c);
% Eliminate negative solutions
check1 = zeros(size(DV_factor));
check_x1 = x1 > check1;
check_x2 = x2 > check1;
x1_mod = x1.*check_x1;
x2_mod = x2.*check_x2;
check_x1_mod = x1_mod >= x2_mod;
check_x2_mod = x1_mod < x2_mod;

dtburn_dt = x1_mod.*check_x1_mod + check_x2_mod.*x2_mod; %Burn fraction for an associated DV_factor
        
% For any burns that can take longer than limit (arent limited by DV2), set
% any values above dtburn_dt_max to dtburn_dt_max 
excess_check = dtburn_dt > dtburn_dt_max;

dtburn_dt = dtburn_dt - excess_check.*dtburn_dt + dtburn_dt_max.*excess_check;



% Adjust DV_factor based on any new dtburn_dt
DV_factor_adj = (polyval(p,dtburn_dt)).^.5;





end
