% Created by Joey Hammond and Colton Hook
% Function coding by Chris Larking, Christian Fuller, Miles Grove, Max McDermott
%
% This code and function sets evaluates a given mission concept and finds
% the optimal mass breakdown and maximum payload operational radius to
% maximize mission success. Success percentile and mass breakdown of
% optimal system is returned

%% Future work

% Coverse about E-prop first, make decision w/ Clarkin (Colton)             Fri 2/26
% Fuel dump functionality                                                   Mon 3/1
% Develop radiator sizing models                                            Mon 3/1

% Ryan:
%   - Find margins and factors for all systems                              Mon 3/1
%   - Add any mass, power draws, and volume additions from SMAP             Mon 3/1
%   - Run sensitivity studies, determine mission concept                    Fri 3/6

% Delageted Work:
%   - Orbits: Expand (and fix) lambert factor model (Joey oversees)         Fri 3/6
%             Expand gravity assist data processing model (Colton oversees) Fri 3/6
%   - Power: Create new power functions/options (Ryan)                      Mon 3/1
%   - Propulsion: Create better structure, volume, etc. (Ryan)              Mon 3/1 
%                 Catolog all options (Ryan)                                Mon 3/1
%                 Get better preposition costs (Ryan)                       Mon 3/1
%                 Get better Xenon density
%   - GNC: Get better funciton of flyby velocity (Ryan)                     Mon 3/1



% Consider coms power draw (distance traveled and what power draw is)








% Dv2  -> Dv1_max and dt1_max
% Dv1 is related to Dv2

% so 
% DV2 can range from 0-DV2
% This coresponds to a DV1_max and a dt1_max
% Knowing engine, DV1 can be reused to find dt


% data shows linear relationship can be assumed, calculate DV1 max and then
% use to find DV1 for any given DV2





%% Inputs



clear
clc
close all

%input orbit data, include numbers and rendevous location for non-instantanous factor
orbitname = 'Orb_1.0x1.15AU_10AULim'; % (From Orbits team)
orbit = [1,1.15]; %AU (From Orbits team)
preposition_DV1 = 3000; %m/s   ---> launch vehicle (From Orbits team)
preposition_DV2 = 300; %m/s    ---> burn 1 (From Orbits team)
LV_mass_capacity = .9; % Fraction of carried launch vehicle mass to capacity (User input)

mass_payload = 600; %Payload mass excluding solar panels and radiators, kg (From SPaM)
flyby_velocity_p = [0,0,10]; % [a,b,c] where flyby velocity (km/s) = ax^2 + bx + c where x is reciprocal of heliocentric range in 1/AU
power_payload = 500; %W (From SPaM)
volume_payload = 10; %m^3 (From SPaM)

R_max = [3, 7]; %Range of heliocentric rendevous design pts, AU (User input)
m_break = [.05,.4]; %Range of mass breakdown (propmass of departure stage/propmass of arrival and departure stages) (User input)
numR2 = 1; % Simulation size (User input)
numMass = 1; % Simulation size (User input)

% prop = [8];
%                   1 - 0 for chemprop (and impulsive) || Thrust for eprop (and non-impulsive), [N]
%                   2 - Thruster dry mass (excluding tank), [kg]
%                   3 - Isp, [s]
%                   4 - Power required, [W]
%                   5 - System Volume (excluding tank), [m^3]
%                   6 - Mixing Ratio (inf for eprop)[O/F]
%                   7 - Oxidizer Density [kg/m^3]
%                   8 - Fuel Density (inf for eprop)[kg/m^3]

XR100 =   [5,  250, 5000, 100000,0,inf,1000,inf]; % XR-100 systems (GUESS IS 1000kg/m^3!!)
XR100_2 = [10, 500, 5000, 200000,0,inf,1000,inf]; %2 XR-100 systems (GUESS IS 1000kg/m^3!!)
% R4D = [0, 3.63, 312, 46]; % 1 R4D system
R4D = [0, 3.63, 312, 0, 0, 1.65, 1440, 880]; % 1 R4D system


% prop_scheme = [preposition_DV2; departure_DV; arrival_DV]
prop_scheme = [R4D;XR100_2;R4D];



%% Current Assumptions (function assumptions not included)

% ----- Preposition -----
% Flight vehicle is Falcon Heavy Expendable
% Flight vehicle uses preposition_DV1
% Mission only viable at end of preposition process
% Max out launch vehicle mass
% No power draw

% -----  Departure  -----
% Solar panel at BOL throughout analysis
% Solar panels unnafected by transfer (AU >= .8)
% Solar panels sized to apogee
% 1 burn

% -----   Arrival   -----
% 1 burn
% No power draw

% -----    Misc     -----
% Standard units: km, yr
% No trajectory adjustment or orbital maintanence burns
% No gravity assists
% Non-instantaneous burns modeled as instantaneous calculated under a safety factor
% Electrical system loss of .775
% Solar panels specific power at 1AU is 120 W/kg


%% Create Propulsion Systems

% Find mass in launch vehicle
[p1,p2] = NonInstantaneousLambert(orbitname);

[m1,V_max] = launchvehicle(preposition_DV1);

% Find mass in preopositioned orbit
preposition_system = prop_scheme(1,:);
R1 = orbit(1);

[mass_array2,power_area2, ~, V2] = prop_sizing1(m1*LV_mass_capacity, 0, R1, preposition_DV2, preposition_system);
m2 = mass_array2(1);


% Simulate a variety of proposed systems
[DV1, DV2, DV1_MAX, DT1, DT2, V, R2,m_break_array,m_test] = propsystemsim(m2, mass_payload, power_payload, prop_scheme, R1, R_max,m_break,numR2,numMass);
V_total = (V + V2 + volume_payload)/V_max;


%% Evaluate Prop System Success Rate
% Determine success rate of each proposed system
maindata = LoadWholeOrbit(orbitname,12); %load in data to matlab once **If iteration used, comment this line, make sure data loaded before while loop**
PosNum = 0; %set to zero means check every orbit position
PercentCoverage = zeros(numR2,numMass); %preallocate coverage matrix



% tic 

for ii = 1:size(DV1,1)
    for jj = 1:size(DV1,2)
        Vsys = V_total(ii,jj);
        DV1sys = DV1(ii,jj); %pull single DV1 for system to check
        DV2sys = DV2(ii,jj); %pull single DV2 for system to check
            DV1_MAXsys = DV1_MAX(ii,jj);
        tburn1 = DT1(ii,jj); %pull single tburn1 for system to check
        tburn2 = DT2(ii,jj); %pull single tburn2 for system to check
        R2sys = R2(ii);      %pull single R2 for system to check (should this be ii or jj?)
        m2 = m_test(ii,jj,1);
        m3 = m_test(ii,jj,2);
        
        tic
        SuccessList = CheckSystem_test(orbitname,PosNum,Vsys,DV1sys,DV2sys,DV1_MAXsys, tburn1,tburn2,R2sys,p1,p2,flyby_velocity_p,maindata,prop_scheme,m2,m3); %produce list of successful ISOs
        PercentCoverage_sys = MixOrbitPositions(SuccessList,SuccessList); %check 2 spacecraft in the same orbit
        [PercentCoverage(ii,jj),~] = BestCombo(PercentCoverage_sys); %checks each spacing for best average coverage, not currently saving spacing value.
        t1 = toc
        
        tic
        SuccessList1 = CheckSystem(orbitname,PosNum,Vsys,DV1sys,DV2sys, tburn1,tburn2,R2sys,p1,p2,flyby_velocity_p,maindata); %produce list of successful ISOs
        PercentCoverage_sys1 = MixOrbitPositions(SuccessList1,SuccessList1); %check 2 spacecraft in the same orbit
        [PercentCoverage1(ii,jj),~] = BestCombo(PercentCoverage_sys1); %checks each spacing for best average coverage, not currently saving spacing value.
        t2 = toc
        
        if ii == 1 && jj == 1
%             fprintf("0%%\n") %indicate loop start
            t_est = toc;
            fprintf("%.1f s Estimated Runtime, %.2f s/system\n",t_est*(numR2*numMass),t_est) %display time taken
            tic
        end
    end
    fprintf("%.1f%%\n",100 * ii/size(DV1,1)) %note progress
end

 t_taken = toc;
 t_taken = t_taken + t_est;
 fprintf("%.1f s Total Runtime Runtime, %.2f s/system\n",t_taken,t_taken/(numR2*numMass)) %display time taken

%% Display Results
results(PercentCoverage,preposition_DV1,preposition_DV2,V_total,DV1,DV2,mass_array2,R2,m_break_array,mass_payload,power_payload,prop_scheme,R1,V_max,V2,volume_payload)


