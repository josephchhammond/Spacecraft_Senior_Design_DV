function [power, mass, area_out] = panel_power(R,area,power)
% This function takes in the distance from the solar panel to the sun, the
% area of the solar panel, and its efficiency and returns the power
% available and mass. Optionally, the power may be passed as the fourth
% input, in which case the function returns a third output of the area
% required to generate that power.
%
% Assumptions:
% - Constant efficiency
% - (1-.775) loss factor (from R.Meisberger's "solarsize.m")
% - Assumes specific power of 120 Watts/kg
% - Assumes rays are all normal to panel (fair for large distance)
%
% Inputs:
% R[1]: Distance from sun to panel, [AU]
% area[1]: Area of solar panel, [m2]
% power[1]: Desired power output of panel, [W] (OPTIONAL)
%
% Outputs:
% power[1]: Power output of panel, [W]
% mass[1]: Mass of panel, [kg]
% area_out[1]: Required panel area, [m2] (OPTIONAL)
%
% Written by: CLarkin
% 
% Edit by Joey -> Took out efficiency from mass calc, specific power should
% have efficiency in it

R = R * 1.496e11;   % Convert from AU to meters

% Constants
P_sun = 3.84e26;    % Total power of sun, [W]
w_kg = 120;         % Specific power of solar array, [W/kg]

% Efficiency loss correction (From R.Meisberger)77
eff = .33 * .775;   % Assume 33% efficiency, 22.5% Loss

switch nargin
    case 2  % User passes area, distance, and efficency to get power out

        % Use inverse-square law to find power at panel
        power = P_sun * area / (4*pi*R^2) * eff; % Array output power, [W]
        area_out = area;        
    case 3  % User passes area, distance, efficiency, and desired power
            % to get required panel area.
        
        % Required panel area
        area_out = power * (4*pi*R^2) / P_sun / eff;  % Panel area, [m2]
end
mass = power / w_kg * (R/1.496e11)^2;  % Mass of array, [kg]
end

