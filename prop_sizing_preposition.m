function [mass_array,P_nom,dt,V] = ...
    prop_sizing_preposition(total_mass, R, dv, prop_system,SMAP)
% This function takes in properties of a preposition system and returns values for
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
% total_mass[1]: Mass of entire stage, [kg]
% power_area[1]: Total area of solar array on payload, [m2]
% R[1]: Distance from sun, [AU]
% dv[1]: Delta-V required for burn, [m/s]
% prop_system[1,4]: Propulsion system characteristics:
%                   1 - 0 for non-impulsive, []|| Thrust for impulsive, [N]
%                   2 - Thruster dry mass, [kg]
%                   3 - Isp, [s]
%                   4 - Power required, [W]
%                   5 - System Volume (m^3)
%                   6 - Mixing Ratio (O/F)
%                   7 - Oxidizer Density (kg/m^3)
%                   8 - Fuel Density (kg/m^3)
% SMAP[3x7]: Allocations external to the propulsion systems ()
%   Rows:
%       1 - Preposition
%       2 - Departure
%       3 - Arrival and payload
%   Columns:
%       1 - Masses [kg]
%       2 - Power required [W]
%       3 - BOL power generation at 1AU [W]
%       4 - BOL radiator wattage dissipation [W]
%       5 - Heat generation by payload [W]
%       6 - Age of spacecraft at current phase [yr]
%       7 - Volume of above components [m^3]

% Outputs:
% mass_array[1,6]:  Mass of different portions of system [kg]
%                   1 - Payload mass (Everything not on this stage), [kg]
%                   2 - Propellant mass, [kg]
%                   3 - Propellant structure mass, [kg]
%                   4 - Power mass, [kg]
%                   5 - Dry mass, [kg]
%                   6 - Total stage mass, [kg]                   
% P_nom[1]: Power generation by solar panels at 1AU added in this stage[W]
% dt[1]: Burn time, [s]
% V[1]: Volume of system, [m^3]
%
% Orignal code written by: CLarkin
% Most edits by : Joey

%% Temporary errors for unadded functionality

%USE UNTIL AGING FUNCTIONALITY ADDED TO RADIATORS AND SP
if norm(SMAP(:,6)) ~=0
    disp("Functionality of age dependent solar panel power generation and")
    disp("radiator heat dissapation not currently supported")
    disp("If you want to proceed assuming beginning of life, press ''y''")
    choice = input('Input choice then press enter.','s');
    if isempty(choice)
        choice = 'y';
    end
    while lower(choice(1)) ~= 'y'
        choice = input("Invalid input, enter ''y'' or ''n''",'s');
    end
end

%USE UNTIL RADIATOR SIZING MODEL ADDED
if SMAP(1,4) < SMAP(1,5) || SMAP(1,4) < SMAP(2,5) || SMAP(3,4) < SMAP(3,5)
    disp("Functionality of automatic radiator sizing not currently supported,")
    disp("and current heat exceeds radiator capabilities in at least one stage")
    disp("If you want to proceed anyways, press ''y''")
    choice = input('Input choice then press enter.','s');
    if isempty(choice)
        choice = 'y';
    end
    while lower(choice(1)) ~= 'y'
        choice = input("Invalid input, enter ''y'' or ''n''",'s');
    end
end

%% Input validation to confirm units are m/s
if dv < 100
    disp("DV should be entered in m/s, it looks like you entered km/s")
    disp("If you entered km/s, press ''y''")
    disp("If you entered this value in m/s and don't want it converted,")
    disp("press ''n''")
    choice = input('Input choice then press enter.','s');
    if isempty(choice)
        choice = 'y';
    end
    while lower(choice(1)) ~= 'y' && lower(choice(1)) ~= 'n'
        choice = input("Invalid input, enter ''y'' or ''n''",'s');
    end
    
    if lower(choice) == 'y'
        dv = dv * 1000;
    elseif lower(choice) == 'n'
    else
        disp("Error in prop_sizing2.m")
        disp("Error with dv unit conversion.")
        pause;
    end
end


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
f_inert = 1/exp(dv/(prop_system(3)*g0));

prop_mass = total_mass*(1-f_inert);                               % Propellant Mass [kg]
struct_mass = prop_mass*struct/(1-struct);                        % Prop struct mass [kg]
mass_array(2) = prop_mass;
mass_array(3) = struct_mass;


P_req = prop_system(4) + SMAP(1,2);
P_onboard = sum(SMAP(:,3));
if P_onboard >= P_req
    mass_array(4) = 0;
    V_panels = 0;
    P_nom = 0;
else
    P_extra = P_req-P_onboard;
    [~,mass_array(4),power_area,V_panels] = panel_power(R, [], P_extra);
    P_nom = P_extra*R^2;
end

thruster_mass = prop_system(2);
dry_mass = thruster_mass + struct_mass + mass_array(4) + SMAP(1,1);         % Dry mass [kg]
mass_array(5) = dry_mass;
mass_array(6) = total_mass;                                             % Intial mass, [kg]

payload_mass = total_mass - dry_mass - prop_mass; % Payload Mass [kg]
mass_array(1) = payload_mass;

if prop_system(1) > 0
    dt = g0*prop_system(3)*prop_mass/prop_system(1); %burn time (s)
else
    dt = 0; %burn time (s)
end

if mass_array(1) < 0
    disp('Error: Infeasible DV for prop system')
    pause;
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

V = V_panels + V_prop + V_sys + SMAP(1,7);

if ~isreal(V) || isnan(V) || V~=norm(V)
    error('Unreal volume of stage')
end

end