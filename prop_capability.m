function [out1] = prop_capability(mass_array, prop_system1, prop_system2)
% This function calculates performance for a given system assuming the
% following:

% 2 stages are present, Eprop and Chem
% Any % of DV1 can be burned
% Any % of DV2 can be burned
% Additional fuel is "dumped" pre-DV1 with an Isp of 300s, with chem dump
% first

% Inputs:
% mass_array[1,6]:  Mass of different portions of system [kg]
%                   1 - Payload mass (Everything not on this stage), [kg]
%                   2 - Eprop Propellant mass, [kg]
%                   3 - Chemprop Propellant mass, [kg]
%                   4 - Eprop Propellant structure mass, [kg]
%                   5 - Chemprop Propellant structure mass, [kg]
%                   6 - Power mass, [kg]
%                   7 - Dry mass, [kg]
%                   8 - Total stage mass, [kg] 
%
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
% Written by: Joseph Hammond


g =9.81; %m/s^2
Isp_dump = 300; %s
Isp1 = prop_system1(3);
Isp2 = prop_system2(3);
T1 = prop_system1(1);

m_total = mass_array(8);
m_prop1 = mass_array(2);
m_prop2 = mass_array(4);
m3 = m_total-m_prop1-m_prop2;%inert mass (mass after fuel is burned)


n = 1;
for ii = 1:100
    for jj = 1:100
        m_burned1 = m_prop1*ii/100;
        m_burned2 = m_prop2*jj/100;
        m_dumped1 = m_prop1 - m_burned1;
        m_dumped2 = m_prop2 - m_burned2;
        
        m0 = m_total;%mass before fuel dump
        m1 = m_total - m_dumped1 - m_dumped2;%mass after fuel dump
        m2 = m1 - m_burned1;%mass after departure
        
        f_dump = m1/m0;
        f_1 = m2/m1;
        f_2 = m3/m2;
        
        DV1_actual(n) = -Isp1*g*log(f_1);
        DV2_actual(n) = -Isp2*g*log(f_2);
        dt_actual(n) = g*Isp1*m_burned1/T1;
        DV_dump(n) = -Isp_dump*g*log(f_dump);
        
        n = n+1;
    end
end

num = 1;
for ii = 1:100
    for jj = 1:100
        
    end
end

figure
title('DV1 vs DV2')
plot(DV1_actual,DV2_actual,'x')
xlabel('DV1 (m/s)')
ylabel('DV2 (m/s)')

figure
title('DVdump vs DV2')
plot(DV_dump,DV2_actual,'x')
xlabel('DV1 (m/s)')
ylabel('DV2 (m/s)')

figure
title('Time vs DV1')
plot(dt_actual,DV1_actual,'x')
xlabel('Time (s)')
ylabel('DV1 (m/s)')

out1 = 1;
end