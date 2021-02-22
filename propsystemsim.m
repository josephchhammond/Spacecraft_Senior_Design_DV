function [DV1, DV2, DT1, DT2, V, R2, m_break] = propsystemsim(m_total,mass_payload, power_payload, prop_scheme, R1, R_max,m_break,numR2,numMass);
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
% Assumptions:
%   Minimum mass breakdown is 10% and maximum mass breakdown is 90%
%
% Outputs:
% Arrays in nxn format with R varied in the 1st dimension and m_breakdown
% varied in the 2nd dimension
%   DV1 = departure burn capability of system, km/s
%   DV2 = arrival burn capability of system, km/s
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
DT1 = zeros(numR2,numMass);
DT2 = zeros(numR2,numMass);
V = zeros(numR2,numMass);
% Determine solar panels required for EP 
EP_power = prop_scheme(2,4);
[~, EP_SP_mass, EP_SP_area,V_EP_panels] = panel_power(R1,[],EP_power);


for ii = 1:numR2
    % Determine mass of solar panels and payload
    R2_temp = R2(ii);
    [~,payload_SP_mass,payload_SP_area,V_payload_panels] = panel_power(R2_temp, [], power_payload);
    m4 = payload_SP_mass+mass_payload;
    
    %Determine variable mass in system
    if payload_SP_mass > EP_SP_mass
        m_available = m2 - m4;
    else
        m_available = m2 - mass_payload - EP_SP_mass;
    end
    
    
    for jj = 1:numMass
        % Determine mass breakdown
        m_break_temp = m_break(jj);
        m3 = m4 + m_break_temp*m_available;
        
        if (m4 > m3) || (m3 > m2)
            error('Mass breakdown error')
        end
        % Create system and develop performance outputs  
        [~,~, dv2,dt2,V2] = prop_sizing2(m4, m3, payload_SP_area, R2_temp, prop_scheme(3,:));
        [~,~, dv1,dt1,V1] = prop_sizing2(m3, m2, payload_SP_area, R1, prop_scheme(2,:));
        
        
        v_total = V1 + V2 + V_payload_panels;
        %Check for any unrealistic answers (negative or imaginary)
        if dv2 ~= norm(dv2) || dv1 ~= norm(dv1) || v_total~= norm(v_total)
            dv2 = 0;
            dt2 = 0;
            dv1 = 0;
            dt1 = 0;
            v_total = 0;
        end
            
        
        DV1(ii,jj) = dv1;
        DV2(ii,jj) = dv2;
        DT1(ii,jj) = dt1;
        DT2(ii,jj) = dt2; 
        V(ii,jj) = v_total;
        
    end
end
end

