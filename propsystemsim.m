function [DV1, DV2, DT, V, R2, m_break, m_array] = propsystemsim(m_total, mass_payload, power_payload, preposition_DV, prop_scheme, R1, SMAP, R_max,m_break,numR2,numMass)
% This function calculates an array of potential propulsion system concepts
% and returns their associated propulsive behaviors
%
% Inputs:
%   m_total = mass of entire vehicle, kg
%   mass_payload = mass of final payload (excluding solar panel mass), kg
%   power_payload = power required by payload, kg
%   prop_scheme = array of propulsion parameters
%           Rows:
%               prepositioning burn
%               departure burn
%               arrival burn
%           Columns:
%               prop = [thrust/type, drymass, ISP, power];
%               thrust/type: = 0 if chemical/instantaneous thrust, else thrust in N
%               drymass: mass of system independent of power and propellant/structure in kg
%               ISP: specific impulse in s
%               power: power required for operation in W
%   R1 = heliocentric range at departure design point, Au
%   R_max = minimum and maximum design points for payload design point, Au
% 
%
% Outputs:
% Arrays in nxn format with R varied in the 1st dimension and m_breakdown
% varied in the 2nd dimension
%   DV1 = departure burn capability of system, km/s
%   DV2 = arrival burn capability of system, km/s
%   DV1_MAX = max carrival departure burn capability if all of arrival system fuel is dumped, km/s
%   DT1 = burn time for departure burn, s
%   DT2 = burn time for arrival burn, s
%
% Array is 1xn
%   R2 = heliocentric range for payload design point, Au
%   m_break = breakdown of available mass to burn 2, %
%
% Written by: Joey Hammond


% Establish arrays and scalar
R2 = linspace(R_max(1),R_max(2),numR2);
m_break = linspace(m_break(1),m_break(2),numMass);
m2 =  m_total;
DV1 = zeros(numR2,numMass);
DV2 = zeros(numR2,numMass);
DT = zeros(numR2,numMass);
m_array = zeros(10,numR2,numMass);
V = zeros(numR2,numMass);
m_stage = SMAP(2,1);


% Determine solar panels required for EP 
EP_power = prop_scheme(1,4);
departure_power = EP_power + SMAP(2,2);
[~, EP_SP_mass, ~,~] = panel_power(R1,[],departure_power);
departure_power = sum(SMAP(:,3)) / R1^2; % Onboard solar panels for departure and on payload
if EP_power>departure_power
    disp("Departure power draw exceeds solar panels");
    pause;
end

% Determine radiators and heat loads present
departure_radiators = sum(SMAP(:,4));
departure_heat = SMAP(2,5) + prop_scheme(1,9);
if departure_heat>departure_radiators
    disp("Departure heat loads exceed radiators onboard, and current functionality doesn't automatically add radiators");
    pause;
end


for ii = 1:numR2
    % Determine mass of solar panels on payload
    R2_temp = R2(ii);
    [~,payload_SP_mass,payload_SP_area,V_payload_panels] = panel_power(R2_temp, [], power_payload);
    
    %Determine variable mass in system
    if payload_SP_mass < EP_SP_mass
        m4 = mass_payload;
        m_available = m2 - m4;
    else
        m4 = payload_SP_mass+mass_payload;
        m_available = m2 - mass_payload - EP_SP_mass;
        disp('Warning - Additional Solar Panels Required in System')
        pause;
    end
    
    
    for jj = 1:numMass
        % Determine mass breakdown
        m_break_temp = m_break(jj);
        m3 = m4 + m_break_temp*m_available;
        
        if (m4 > m3) || (m3 > m2)
            error('Mass breakdown error')
        end
        
        
        [mass_array,power_area,dv1, dv2, dt,v] = ...
            prop_sizing4(m4, m2, m_stage, m_break_temp, departure_power, R1, prop_scheme(1,:), prop_scheme(2,:),preposition_DV);
        
         
%         if ~isreal(mass_array) || isnan(mass_array)
%             error('Unreal mass array')
        if ~isreal(power_area) || isnan(power_area)
            error('Unreal power area')
        elseif ~isreal(dv1) || isnan(dv1)
            error('Unreal DV1')
        elseif ~isreal(dv2) || isnan(dv2)
            error('Unreal DV2')
        elseif ~isreal(dt) || isnan(dt)
            error('Unreal burn time')
        end
        
        DV1(ii,jj) = dv1;
        DV2(ii,jj) = dv2;
        DT(ii,jj) = dt;
        V(ii,jj) = v;
        m_array(:,ii,jj) = [mass_array];
        
    end
end
end

