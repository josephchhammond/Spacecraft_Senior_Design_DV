function [DV_adj] = DV_adjustment(DV,p,dtburn,dttransfer)
% This function adjusts DV based on non-instantaneous behavior
% 
% Inputs:
% DV - Instantaneous Delta V in need of adjustment
% p - Coefficients for adjustment equation for noninstantaneous behavior
% dtburn_dt - Fraction of transfer spent in burn
% 
% Output:
% DV_adj - Noninstantaneous Delta V

%DV_adj = DV * (polyval(p,dtburn/dttransfer)).^.5; %original


DV_adj = DV .* (polyval(p,dtburn./dttransfer)).^.5; %changed to this by colton on 2/8/21
end
