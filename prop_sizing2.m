function [mass_array,power_area, dv,dt,V] = ...
    prop_sizing2(payload_mass, m0, power_area, R, prop_system)
% This function takes in properties of the system and returns values for
% the mass of the propulsion system, as well as the necessary total solar
% array mass.
%
% Assumptions:
% - Impulsive burn
% - "Payload mass" is the mass of everything not included in the burn, i.e.
% spacecraft as well as other propulsion systems.
% - Constant thrust, no transients on startup/shutdown
%   
% Inputs:
% payload_mass[1]: Mass of payload, [kg]
% m0[1]: Total initial mass of stage, [kg]
% power_area[1]: Total area of solar array on payload, [m2]
% R[1]: Distance from sun, [AU]
% m0[1]: Initial mass of spacecraft, [kg]
% prop_system[1,4]: Propulsion system characteristics:
%                   1 - 0 for non-impulsive, []|| Thrust for impulsive, [N]
%                   2 - Thruster dry mass, [kg]
%                   3 - Isp, [s]
%                   4 - Power required, [W]
%                   5 - System Volume (m^3)
%                   6 - Mixing Ratio (O/F)
%                   7 - Oxidizer Density (kg/m^3)
%                   8 - Fuel Density (kg/m^3)
%
% Outputs:
% mass_array[1,6]:  Mass of different portions of system [kg]
%                   1 - Payload mass (Everything not on this stage), [kg]
%                   2 - Propellant mass, [kg]
%                   3 - Propellant structure mass, [kg]
%                   4 - Power mass, [kg]
%                   5 - Dry mass, [kg]
%                   6 - Total stage mass, [kg]                   
% power_area[1]: Total area of solar array required, [m2]
% dt[1]: Burn time, [s]
% V[1]: Volume of system, [m^3]
%
% Written by: CLarkin

% Constants
g0 = 9.81;      % Gravitational acceleration, [m/s2]


%% Input validation to confirm instantaneous-compatible prop system
if prop_system(1) > 0       % E-Prop thruster
    % Input validation for Chemical Isp
    if prop_system(3) > 500 && prop_system(3) <= 5000
        % Pass
    elseif prop_system(3) > 5000
        disp("Error in prop_sizing1.m")
        disp("Invalid Isp for chemical (too high)")
        pause;
    else
        disp("Error in prop_sizing1.m")
        disp("Invalid Isp for chemical (too low)")
        pause;
    end
    struct = .05;   % Ratio of of tank mass for added prop and tank mass [unitless]
    dm = prop_system(1) / (g0*prop_system(3));


elseif prop_system(1) == 0   % Chem-Prop thruster
    % Input validation for Chemical Isp
    if prop_system(3) > 200 && prop_system(3) <= 500
        % Pass
    elseif prop_system(3) > 500
        disp("Error in prop_sizing1.m")
        disp("Invalid Isp for e-prop (too high)")
        pause;
    else
        disp("Error in prop_sizing1.m")
        disp("Invalid Isp for e-prop (too low)")
        pause;
    end
    struct = .15;   % Ratio of of tank mass for added prop and tank mass [unitless]
    dm = inf;

    
elseif prop_system(1) < 0   % Catch for anomaly
    disp("Error in prop_sizing1.m")
    disp("Invalid thrust level (Negative)")
    pause;
else
end


%% Mass assignments
mass_array(1) = payload_mass;                           % Payload Mass [kg]
mass_array(6) = m0;                                     % Initial mass [kg]
    
[~,power_mass1,~,V_panels1] = panel_power(R, power_area);                         % Power mass on vehicle [kg]
[~,power_mass2,power_area2,V_panels2] = panel_power(R, [], prop_system(4));        % Power mass required at stage [kg]
if power_mass2 >= power_mass1
    mass_array(4) =  power_mass2 - power_mass1;
    power_area = power_area2;
    V_panels = V_panels2-V_panels1;
else
    mass_array(4) = 0;
    V_panels = 0;
end


mass_array(2) = (1-struct) * (m0 - (mass_array(1) + prop_system(2) + mass_array(4))); % Propellant mass [kg]
mass_array(3) = (struct) * (m0 - (mass_array(1) + prop_system(2) + mass_array(4)));   % Propellant mass [kg]
mass_array(5) = prop_system(2) + mass_array(3) + mass_array(4);                % Dry mass [kg]

% Delta-V calc
dv = prop_system(3)*g0*log(mass_array(6)/(mass_array(6)-mass_array(2)));% Delta-V [km/s]

if prop_system(1) > 0
    dt = g0*prop_system(3)*mass_array(2)/prop_system(1); %burn time (s)
else
    dt = 0; %burn time (s)
end

%% Volume Assignments

O_F = prop_system(6); % Mixing Ratio (O/F)
rho_O = prop_system(7); % Oxidizer Density (kg/m^3)
rho_F = prop_system(8); % Fuel Density (kg/m^3)

m_O = mass_array(2)*(1-1/(O_F+1)); % Oxidizer mass (kg)
m_F = mass_array(2)/(O_F+1); % Fuel mass (kg)


V_O = m_O/rho_O; % Oxidizer volume (m^3)
V_F = m_F/rho_F; % Fuel volume (m^3)
V_sys = prop_system(5);
V_prop = V_O+V_F;

V = V_panels + V_prop + V_sys;

if ~isreal(V) || isnan(V) || V~=norm(V)
    error('Unreal volume of stage')
end

end