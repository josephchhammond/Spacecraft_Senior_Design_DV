% Created by Joey Hammond and Colton Hook
% Function coding by Chris Larking, Christian Fuller, Miles Grove, Max McDermott
%
% This code and function sets evaluates a given mission concept and finds
% the optimal mass breakdown and maximum payload operational radius to
% maximize mission success. Success percentile and mass breakdown of
% optimal system is returned
%% Inputs


clear
clc
close all

%input orbit data, include numbers and rendevous location for non-instantanous factor
orbitname = 'Orb_1.0x1.15AU_10AULim';
orbit = [1,1.35]; %AU
preposition_DV1 = 3000; %m/s   ---> launch vehicle
preposition_DV2 = 300; %m/s   ---> burn 1

mass_payload = 600; %kg
flyby_velocity_p = [0,0,10]; % [a,b,c] where flyby velocity (km/s) = a^2x + bx + c where x is reciprocal of heliocentric range in 1/AU
power_payload = 500; %W
R_max = [3, 11]; %Range of heliocentric rendevous design pts, AU
m_break = [.04,.6]; %Range of mass breakdown (propmass of departure stage/propmass of arrival and departure stages)

% prop = [thrust/type, drymass, ISP, power];
% thrust/type: = 0 if chemical/instantaneous thrust, else thrust in N
% drymass: mass of system independent of power and propellant/structure in kg
% ISP: specific impulse in s
% power: power required for operation in W

XR100 = [5, 250, 5000, 100000]; %2 XR-100 systems
XR100_2 = [10, 500, 5000, 200000]; %2 XR-100 systems
% R4D = [0, 3.63, 312, 46]; % 1 R4D system
R4D = [0, 3.63, 312, 0]; % 1 R4D system


% prop_scheme = [preposition_DV2, departure_DV, arrival_DV] ?Can we remove?
prop_scheme = [R4D;XR100_2;R4D];

% Size of simulation
numR2 = 9;
numMass = 9;

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
%% Create Propulsion Systems

%preliminary code for itteration
% improvement = 1;
% maxval_prev = 0;
% TOL = 1e-2;
% count = 1;
% maindata = LoadWholeOrbit(orbitname,12); %load in data to matlab once
%while improvement>TOL && count<10

% Find mass in launch vehicle
[p1,p2] = NonInstantaneousLambert(orbitname);
[m1] = launchvehicle(preposition_DV1);

% Find mass in preopositioned orbit
preposition_system = prop_scheme(1,:);
R1 = orbit(1);

[mass_array2,power_area2, ~] = prop_sizing1(m1, 0, R1, preposition_DV2, preposition_system);
m2 = mass_array2(1);


% Simulate a variety of proposed systems
[DV1, DV2, DT1, DT2, R2,m_break_array] = propsystemsim(m2, mass_payload, power_payload, prop_scheme, R1, R_max,m_break,numR2,numMass);

% n = 1;
% for ii = 1:size(DV1,1)
%     for jj = 1:size(DV1,2)
%         DV1_plot(n) = DV1(ii,jj)/1000;
%         DV2_plot(n) = DV2(ii,jj)/1000;
%         n = n+1;
%     end
% end
% plot(DV1_plot,DV2_plot,'x')
% xlabel('Burn 1 (km/s)')
% ylabel('Burn 2 (km/s)')
% title('Preprocessed DV Capabilities (Blind to dt)')


%% Evaluate Prop System Success Rate
% Determine success rate of each proposed system
maindata = LoadWholeOrbit(orbitname,12); %load in data to matlab once **If iteration used, comment this line, make sure data loaded before while loop**
PosNum = 0; %set to zero means check every orbit position
PercentCoverage = zeros(numR2,numMass); %preallocate coverage matrix
tic 

fprintf("0%%\n") %indicate loop start
for ii = 1:size(DV1,1)
    for jj = 1:size(DV1,2)
        DV1sys = DV1(ii,jj); %pull single DV1 for system to check
        DV2sys = DV2(ii,jj); %pull single DV2 for system to check
        tburn1 = DT1(ii,jj); %pull single tburn1 for system to check
        tburn2 = DT2(ii,jj); %pull single tburn2 for system to check
        R2sys = R2(ii);      %pull single R2 for system to check (should this be ii or jj?)
        SuccessList = CheckSystem(orbitname,PosNum,DV1sys,DV2sys,tburn1,tburn2,R2sys,p1,p2,flyby_velocity_p,maindata); %produce list of successful ISOs
        PercentCoverage_sys = MixOrbitPositions(SuccessList,SuccessList); %check 2 spacecraft in the same orbit
%         PercentCoverage(ii,jj) = mean(diag(PercentCoverage_sys)); %assume only evaluating 1 spacecraft (average over whole orbit)
        %PercentCoverage_sys = [PercentCoverage_sys(:,7:12),PercentCoverage_sys(:,1:6)];
        %PercentCoverage(ii,jj) = mean(diag(PercentCoverage_sys)); %assume evaluating 2 spacecraft in 180deg displacement(average over whole orbit)
        [PercentCoverage(ii,jj),~] = BestCombo(PercentCoverage_sys); %checks each spacing for best average coverage, not currently saving spacing value.
    end
    fprintf("%.1f%%\n",100 * ii/size(DV1,1)) %note progress
end
t_taken = toc;
fprintf("%.1f s taken, %.2f s/system\n",t_taken,t_taken/(numR2*numMass)) %display time taken


%prelimnary iteration code
% fprintf("Iteration %i, best coverage: %.2f%%\n",count,maxval)
% R_max = [R2(max(ii-2,1)),R2(min(ii+2,length(R2)))];
% m_break = [m_break_array(max(jj-2,1)),m_break_array(min(jj+2,length(m_break_array)))];
% improvement = abs(maxval_prev-maxval)/maxval;
% maxval_prev = maxval;
% count = count+1;
%pause
%end
%%

results(PercentCoverage,preposition_DV1,preposition_DV2,DV1,DV2,mass_array2,R2,m_break_array,mass_payload,power_payload,prop_scheme,R1)


