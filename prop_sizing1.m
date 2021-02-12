function [mass_array,power_area,dt] = ...
    prop_sizing1(total_mass, power_area, R, dv, prop_system)
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
% power_area[1]: Total area of solar array on payload, [m2]
% R[1]: Distance from sun, [AU]
% dv[1]: Delta-V required for burn, [m/s]
% prop_system[1,4]: Propulsion system characteristics:
%                   1 - 0 for non-impulsive, []|| Thrust for impulsive, [N]
%                   2 - Thruster dry mass, [kg]
%                   3 - Isp, [s]
%                   4 - Power required, [W]
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
%
% Written by: CLarkin
% Past edits: Joey

% Input validation to confirm units are m/s
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


% Input validation to confirm instantaneous-compatible prop system
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



% Mass assignments
f_inert = 1/exp(dv/(prop_system(3)*g0));
mass_array(2) = total_mass*(1-f_inert);                                 % Propellant Mass [kg]
mass_array(3) = mass_array(2)*struct/(1-struct);                        % Prop struct mass [kg]

[~,power_mass1,~] = panel_power(R, power_area);                         % Power mass on vehicle [kg]
[~,power_mass2,power_area2] = panel_power(R, [], prop_system(4));       % Power mass required at stage [kg]
if power_mass2 >= power_mass1
    mass_array(4) =  power_mass2 - power_mass1;
    power_area = power_area2;
else
    mass_array(4) = 0;
end

mass_array(5) = prop_system(2) + mass_array(3) + mass_array(4);         % Dry mass [kg]
mass_array(6) = total_mass;                                             % Intial mass, [kg]

mass_array(1) = mass_array(6) - mass_array(5) - mass_array(2); % Payload Mass [kg]

if prop_system(1) > 0
    dt = g0*prop_system(3)*mass_array(2)/prop_system(1); %burn time (s)
else
    dt = 0; %burn time (s)
end

if mass_array(1) < 0
    disp('Error: Infeasible DV for prop system')
    pause;
end



% % Feasibility Condition check
% f_inert = payload_mass/mass_array(6);
% if prop_system(3) <= dv/(log(1/f_inert)*g0)
%     disp("Error in chem_prop_sizing.m")
%     disp("Error: Non-feasible condition.")
%     pause;
% else    % If it passes feasibility condition, run values
%     mass_array(2) = mass_array(6) - mass_array(5);
%     
%     % Check if power is sufficient
%     if panel_power(R, power_area) >= prop_system(4)
%         % Power is sufficient, get array mass.
%         disp("Current solar array size is sufficient for prop system.")
%         [~,mass_array(4)] = panel_power(R, power_area);
%     elseif panel_power(R, power_area) < prop_system(4)
%         % Power isn't sufficient, resize and get array mass.
%         disp("Current solar array size is not sufficient for prop system.")
%         disp("Upsizing solar array.")
%         [~,mass_array(4), power_area] = panel_power(R, power_area, ...
%         prop_system(4));
%     else
%         % Anomaly catch
%         disp("Error in chem_prop_sizing.m")
%         disp("There was a problem calculating solar array area.")
%         pause;
%     end
% end

% disp(' ')
% fprintf('Payload mass: %f kg\n', mass_array(1))
% fprintf('Propellant mass: %f kg\n', mass_array(2))
% fprintf('Propellant structure mass: %f kg\n', mass_array(3))
% fprintf('Power mass: %f kg\n', mass_array(4))
% fprintf('Dry mass: %f kg\n', mass_array(5))
% fprintf('Total stage mass: %f kg\n', mass_array(6))
% fprintf('Total solar array area required: %f m2\n', power_area)
% fprintf('Total burn time is: %f m2\n', dt)
% disp('~~~~~')
end