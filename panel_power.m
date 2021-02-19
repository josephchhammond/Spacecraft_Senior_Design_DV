function [power, mass, area_out, volume] = panel_power(R,area,power)
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




% Constants
% AU = 1.496e11;              % Convert from AU to meters
% P_sun = 3.84e26;            % Total power of sun, [W]
% Efficiency loss correction (From R.Meisberger)
eff = .775;                 % Assume 22.5% loss in system
w_kg_1AU = 120;             % Specific power of solar array at 1AU, [W/kg]
w_density_1AU = 20000/7;    % [W/m^3] Volume of 1 20kW array is 0.6096x6m cylinder = 7m^2
w_A_1AU = 20000/(6*13.7);   %[W/m^2]

%Add a cutoff for minimum radius, assumption is if we are under 1AU we can
%protect solar panels by angleing them so incident radiation is that of 1AU
if R < 1
    R = 1;
end

%P_1AU = P_sun / (4*pi*AU^2); % Solar irradiance at 1AU
w_kg = w_kg_1AU / R^2 * eff;
w_density = w_density_1AU / R^2 * eff;
w_A = w_A_1AU / R^2 * eff;
switch nargin
    case 2  % User passes area and distance to get power out

        % Use inverse-square law to find power at panel
        area_out = area;
        power = w_A * area * eff;
        
                % power = P_1AU*eff/R^2; % Array output power, [W]

    case 3  % User passes area, distance, efficiency, and desired power
            % to get required panel area.
        
        % Required panel area
        area_out = power/(eff*w_A);  % Panel area, [m2]
end

mass = power / w_kg;        % Mass of array, [kg]
volume = power / w_density; % Volume of stored array, [kg]
end

