function [mass_array,power_area,dv1, dv2, dt,V] = ...
    prop_sizing3(payload_mass, m0, m_break, power_area, R, prop_system1, prop_system2)
% This function calculates performance for a given departure and arrival
% system on a thrust plate.
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
% m_break[1]: Breakdown of propellant mass (fraction on departure system) []
% power_area[1]: Total area of solar array on payload, [m2]
% R[1]: Distance from sun, [AU]
% m0[1]: Initial mass of spacecraft, [kg]
% prop_system1[1,8]: Departure eprop propulsion system characteristics:
% prop_system2[1,8]: Arrival chempop propulsion system characteristics:
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
%                   2 - Eprop Propellant mass, [kg]
%                   3 - Chemprop Propellant mass, [kg]
%                   4 - Eprop Propellant structure mass, [kg]
%                   5 - Chemprop Propellant structure mass, [kg]
%                   6 - Power mass, [kg]
%                   7 - Dry mass, [kg]
%                   8 - Total stage mass, [kg] 

% power_area[1]: Total area of solar array required, [m2]
% dv1[1]: Departure burn DV (assumed both full chemprop and eprop tanks), [m/s]
% dv2[1]: Arrival burn DV (assumed  full chemprop tank and empty eprop tank), [m/s]
% dt[1]: Burn time, [s]
% V[1]: Volume of system, [m^3]
%
% Written by: Joseph Hammond


%% Input validation to confirm instantaneous-compatible prop system
if prop_system1(1) > 0       % E-Prop thruster
    % Input validation for Chemical Isp
    if prop_system1(3) > 500 && prop_system1(3) <= 5000
        % Pass
    elseif prop_system1(3) > 5000
        error("Invalid Isp for E-prop (too high)")
        pause;
    else
        error("Invalid Isp for E-prop (too low)")
        pause;
    end
    struct1 = .05;   % Ratio of of tank mass for added prop and tank mass [unitless]
    dm1 = prop_system1(1) / (g0*prop_system1(3));
else
    error('Departure burn must be E-prop')
end



if prop_system2(1) == 0   % Chem-Prop thruster
    % Input validation for Chemical Isp
    if prop_system2(3) > 200 && prop_system2(3) <= 500
        % Pass
    elseif prop_system2(3) > 500
        error("Invalid Isp for chem-prop (too high)")
    else
        error("Invalid Isp for chem-prop (too low)")
    end
    struct2 = .15;   % Ratio of of tank mass for added prop and tank mass [unitless]
    dm2 = inf;
else
    error('Departure burn must be Chem-prop')
end


%% Mass assignments
mass_array(1) = payload_mass;                           % Payload Mass [kg]
mass_array(8) = m0;                                     % Initial mass [kg]

P_prop = prop_system1(4) + prop_system2(4);
[~,power_mass1,~,V_panels1] = panel_power(R, power_area);                         % Power mass on vehicle [kg]
[~,power_mass2,power_area2,V_panels2] = panel_power(R, [], P_prop);        % Power mass required at stage [kg]
if power_mass2 >= power_mass1
    mass_panels =  power_mass2 - power_mass1;
    power_area = power_area2;
    V_panels = V_panels2-V_panels1;
else
    mass_panels = 0;
    V_panels = 0;
end
mass_array(6) = mass_panels;

thruster_mass = prop_system1(2) + prop_system2(2); % Mass of both thrusters/thruster plate [kg]
m_available = m0 - payload_mass - mass_panels - thruster_mass; % Propellant and tank mass available for both burns [kg]

m_available1 = m_available * m_break; % Propellant and tank mass of departure burn [kg]
m_available2 = m_available - m_available1; % Propellant and tank mass of arrival burn [kg]

m_prop1 = (1-struct1) * m_available1; % Eprop Propellant mass, [kg]
m_tank1 = m_available1 - m_prop1; % Eprop Propellant structure mass, [kg]
m_prop2 = (1-struct2) * m_available2; % Chemprop Propellant mass, [kg]
m_tank2 = m_available2 - m_prop2; % Chemprop Propellant structure mass, [kg]
mass_array(2) = m_prop1;
mass_array(4) = m_tank1;
mass_array(3) = m_prop2;
mass_array(5) = m_tank2;

m_dry = thruster_mass + mass_panels + m_tank1 + m_tank2; % Stage Dry Mass, [kg]
mass_array(7) = m_dry;


%% Maneuvers
g0 = 9.81;      % Gravitational acceleration, [m/s2]
Isp1 = prop_system1(3); %Specific impulse of departure [s]
T1 = prop_system1(1); %Thrust of departure [N]
f_inert1 = (payload_mass - m_prop1)/payload_mass; % Inert mass fraction of departure, [m/s2]

dv1 = Isp1*g0*log(f_inert1);
dt = g0*Isp1*m_prop1/T1;

Isp2 = prop_system2(3); %Specific impulse of departure [s]
f_inert2 = (payload_mass-m_prop1-m_prop2)/(payload_mass-m_prop1); % Inert mass fraction of arrival, [m/s2]

dv2 = Isp2*g0*log(f_inert2);
%% Volume Assignments

% Departure burn
O_F = prop_system1(6); % Mixing Ratio (O/F)
rho_O = prop_system1(7); % Oxidizer Density (kg/m^3)
rho_F = prop_system1(8); % Fuel Density (kg/m^3)

m_O = m_prop1*(1-1/(O_F+1)); % Oxidizer mass (kg)
m_F = m_prop1/(O_F+1); % Fuel mass (kg)


V_O = m_O/rho_O; % Oxidizer volume (m^3)
V_F = m_F/rho_F; % Fuel volume (m^3)
V_sys = prop_system1(5);
V_prop1 = V_O+V_F;


% Arrival burn
O_F = prop_system2(6); % Mixing Ratio (O/F)
rho_O = prop_system2(7); % Oxidizer Density (kg/m^3)
rho_F = prop_system2(8); % Fuel Density (kg/m^3)

m_O = m_prop2*(1-1/(O_F+1)); % Oxidizer mass (kg)
m_F = m_prop2/(O_F+1); % Fuel mass (kg)


V_O = m_O/rho_O; % Oxidizer volume (m^3)
V_F = m_F/rho_F; % Fuel volume (m^3)
V_sys = prop_system2(5);
V_prop2 = V_O+V_F;


V = V_panels + V_prop1 + V_prop2 + V_sys;

if ~isreal(V) || isnan(V) || V~=norm(V)
    error('Unreal volume of stage')
end

end