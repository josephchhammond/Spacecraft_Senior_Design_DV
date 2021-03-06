function [Success] = CheckSystem_test(OrbitName,PosNum,Vsys,DV1sys,DV2sys,DV1_MAXsys,tburn1,tburn2,R2sys,p1_adjust,p2_adjust,p_flyby,MainData,prop_scheme,m1,m2)
%Takes the info about which orbits data to check along with prop system values and
%computes a percent success by calling the other functions defined below.
%
%Inputs:
%OrbitName - string, Name of the orbit (the folder that contains the data)
%PosNum    - integer from 0 to 12, Position in the orbit to check, 0 means check every
%               system and return a 7500x12. other numbers mean check specific position only 
%DV1sys    - scalar value in km/s of the system's DV1
%DV2sys    - scalar value in km/s of the system's DV2
%DV1_MAXsys    - scalar value in km/s of the system's DV1 if all of DV2 fuel is dumped
%tburn1    - scalar value in seconds of the burn time associated with DV1
%tburn2    - scalar value in seconds of the burn time associated with DV2
%R2sys     - scalar value in AU of the R2 the system was designed for
%p1_adjsut- 3x1 array of polynoimial coefficients for correcting DV1
%p2_adjsut- 3x1 array of polynoimial coefficients for correcting DV2
%p_flyby- 3x1 array of polynoimial coefficients for finding max flyby speed
%prop_scheme-propulsive data for mission
%
%Outputs:
%Successes - 7500x12 or 7500x1 depending on input "PosNum". each row is an ISO while each
%               collumn is the orbit position
%
%Written by Colton Hook on 2/8/2021 - This is V4


    %unit adjustments
    AU = 1.496e8; %AU to km
    R2sys = R2sys*AU; % adjust to AU
    DV1sys = DV1sys/1000; %convert to km/s
    DV2sys = DV2sys/1000; %convert to km/s
    Success = []; %create empty variable

    
    
    F1 = prop_scheme(2,1); %Departure: F=thrust [N], Isp=specific impulse [s]
    Isp1 = prop_scheme(2,3); %Departure: F=thrust [N], Isp=specific impulse [s]
    Isp2 = prop_scheme(3,3); %Arrival: Isp=specific impulse [s]
    g = 9.81; %m/s
    
    m1_burned = tburn1*F1 / (g*Isp1); %kg
    m2_burned = m2 * (1-exp(-DV2sys*1000/(Isp2*g))); %kg
      
%if input is zero, loop through positions 1 through 12 and save data after each run
%otherwise just check orbit listed
    if PosNum == 0 %set limits to all include all orbit positions
        a = 1;
        b = 12;
    else %set limits to only include the one orbit position
        a = PosNum;
        b = PosNum;
    end
    
for jj = a:b %for each orbit position
    PosNum = jj; %reusing this variable
%check a specific position for all ISOs add it to the output list
    OrbitData = MainData{1,PosNum}; %load in the data for the orbit and the position
    n = length(OrbitData.dv1); %how many ISOs are in this data set?
    ISOsuccess = zeros(n,1); %initalize the ISOsuccess vector as all zeros
    for ii = 1:n %for each ISO
        [DV1data,DV2data,dtdata,R2data] = readISOdata(OrbitData,ii); %pulls the data for the ii ISO
        [pass_fail] = checktransfer(DV1data,DV2data,dtdata,R2data,Vsys,DV1sys,DV2sys,DV1_MAXsys,tburn1,tburn2,R2sys,p1_adjust,p2_adjust,p_flyby,m1,m2,m1_burned,m2_burned,Isp1,Isp2,F1); %checks the ISO data agianst system data
        num_good_transfers = sum(pass_fail,'all');
        if num_good_transfers > 0 %if any transfer to the ISO succeeded, mark that ISO as a success
            ISOsuccess(ii) = 1;
        end
    end
    Success = [Success,ISOsuccess];
end

end

%***This function is no longer used. Data is loaded once and then read repeatedly***
% function [OrbitData] = LoadOrbit(PosNum,OrbitMain)
% %using the name of the orbit and the position number of interest, outputs the relevant
% %orbit data for all the ISOs from that positon on that orbit. for use with "readISOdata"
% %Inputs:
% %PosNum    - integer from 1 to 12, Position in the orbit to check
% %
% %Outputs:
% %OrbitData - the struct containing all of the relevant data from this orbital position
% %
% %Written by Colton Hook on 2/8/2021 - This is V4
% 
%     filestring = "Results_" + PosNum +".mat";
%     OrbitData = load(filestring);
%     
% end

function [DV1data,DV2data,dtdata,R2data] = readISOdata(OrbitData,ISOnum)
%Reads data for 1 ISO for use with "checktransfer"
%
%Inputs:
%OrbitData - struct obtained from "LoadOrbit"
%ISOnum    - number of the ISO to pull values for
%
%Outputs: (all for direct use with checktransfer
%DV1data - 30x30 of DV1s
%DV2data - 30x30 of DV2s
%dtdata  - 1x30 of transfer times
%R2data  - 30x30 of R2s for each transfer
%
%Written by Colton Hook on 2/8/2021 - This is V4

    DV1data = OrbitData.dv1{ISOnum,1};
    DV2data = OrbitData.dv2{ISOnum,1};
    dtdata  = OrbitData.t_transfer{ISOnum,1};
    R2data  = OrbitData.R2{ISOnum,1};
end

function [pass_fail] = checktransfer(DV1data,DV2data,dtdata,R2data,Vsys,DV1sys,DV2sys,DV1_MAXsys,tburn1,tburn2,R2sys,p1_adjust,p2_adjust,p_flyby,m1,m2,m1_burned,m2_burned,Isp1,Isp2,F1)
% Converts the provided orbits data to a pass fail metric based on the system parameters
% Calculation works backwards in burns to evaluate feasibility. General
% steps of model are as follows:
%
%   1. Check if our system can make our arrival burn. One of three
%   possibilities occurs: 
%       a) We cant make the arrival burn even at a perfect
%       instantaneous departure
%       b) We can make the arrival under a limited amount of
%       "non-instantaenty", and we set an appropriate burn time limit
%       c) We can make the arrival past burn time limit, so we dump
%       extra fuel before departure
%
%   2. For a given burn time, we find how much DV we get out of our
%   departure system and how much fuel has to be dumped
%
%   3. For both fuel dumps, we find how much extra DV we get before we
%   depart
%
%   4. We check if our system can make the arrival and departure
%   requirements for our given burns and dumps
%
%   5. Check if volume and power limitations are met
%
% Inputs:
% DV1data - 30x30 of the DV1s for the ISO
% DV2data - 30x30 of the DV2s for the ISO
% dtdata  -  1x30 of the dts for the ISO, each collumn of the 30x30s have same dt
% R2data  - 30x30 of the R2 for the ISO
% DV1sys  - scalar value of the system's DV1
% DV2sys  - scalar value of the system's DV2
% tburn1  - scalar value of the burn time associated with DV1
% tburn2  - scalar value of the burn time associated with DV2
% R2sys   - scalar value of the R2 the system was designed for
% p1_adjsut- 3x1 array of polynoimial coefficients for correcting DV1
% p2_adjsut- 3x1 array of polynoimial coefficients for correcting DV2
% p_flyby- 3x1 array of polynoimial coefficients for finding max flyby speed
%
% Outputs:
% pass_fail 30x30 logical array 1 means success, 0 means failure for that specific
% transfer, if a single entry is 1, then the ISO can be reached with the system
%
%Written by Colton Hook  and Joey Hammond on 2-3/2021 - This is V6





%  %Inputs to fix and input
%     g = 9.81; %m/s^2
%     Isp_dump1 = 300; %s
%     Isp_dump2 = 300; %s
%     T1
%     Isp1
%     Isp2
%     m_prop1;
%     m_prop2;
%     m3 ;%inert mass (mass after fuel is burned)
%     DV2sys;
%     dtburn_dt_max = 0.5; %change to higher fidelity model later *ASSUMPTION ALERT*
    
 

% Convert dtdata into a 30x30 by extending it down each column
    dtdata_ = repmat(dtdata,30,1);

   
 
% DV2 (fourth burn) 
    % Reduce arrival requirements based on flyby velocity
    maxflyby = polyval(p_flyby,R2data);
    DV2adj_ = max(DV2data - maxflyby,0);
    
    % Find how many times we can cover the arrival requirements
    DV2_factor = DV2sys./DV2adj_;
    
    % Find how much burntime we can allow for this system and actual DV
    % requirements for our arrival stage
        % bad systems have dtburn_dt < 0
        % good systems should have dtburn_dt 0-.5
        % really good systems have dtburn_dt .5 
    [dtburn_dt,DV2_factor_adj] = DV_adjustment2(DV2_factor,p,dtburn_dt_max);    % Iterate fuel dumping and burn time calcs until convergence
    DV2_req = DV2_factor_adj.*DV2adj_; %Arrival DV specific to this transfer
    dtburn = dtdata_.*dtburn_dt;
    % Find fuel dumped
    f_2 = exp(DV2_req/(g*-Isp2));
    m_burned2 = f_2*m3/(1-f_2);
    m_dump2 = m_prop2 - m_burned2;

% DV1 (third burn)
    % Based on allowable burn time, find burned fuel on departure stage and
    % how much is saved
    m_burned1 = dtburn*T1/(g*Isp1);
    m_dump1 = m_prop1 - m_burned1;
    % Find how much departure DV we get, adjusted for non-instantaety
    f_1 = (m3 + m_burned2)/(m3 + m_burned1 + m_burned2);
    DV1_noninst = -Isp1*g*log(f_1);
    DV1_factor = DV_adjustment(1,p1_adjust,dtburn_dt,1);
    DV1_noninst_adj = DV1_noninst ./ DV1_factor;
    
% DV2 dump (second burn)
    % Find how much DV we get out of dumping extra arrival fuel before
    % departure
    m_post_dump2 = m3 + m_burned1 + m_burned2;
    m_pre_dump2 = m_post_dump2 + m_dump2;
    f_dump2 = m_post_dump2/m_pre_dump2;
    DV_dump2 = -Isp_dump2*g*log(f_dump2);
    
% DV1 dump (first burn)
    % Find how much DV we get out of dumping extra departure fuel before
    % departure
    m_post_dump1 = m_pre_dump2;
    m_pre_dump1 = m_post_dump1 + m_dump1;
    f_dump1 = m_post_dump1/m_pre_dump1;
    DV_dump1 = -Isp_dump1*g*log(f_dump1);
    
% Tally total departure DV output
    DV1_total = DV_dump1 + DV_dump2 + DV1_noninst_adj;
    
    
    %Criteria to meet
    crit1 = R2data<= R2sys; %we are capable than going further than this ISO
    crit2 = DV2_factor >= 1; %Check if DV2 is met
    crit3 = DV1_total >= DV1data; %Check if DV1 is met under burntime limitations
    crit4 = Vsys <= 1;
    
    pass_fail = crit1 & crit2 & crit3 & crit4;
end

%---------------------------Change Log---------------------------
%V1 - created 2/2/2021
%V2 - updated 2/3/2021
%   updated crit2 and crit3 in "checktransfer" (they were backwards)
%   added conversion of R2sys to km in "CheckSystem"
%   clarified units for inputs to "CheckSystem"
%V3 - updated 2/5/2021
%   updated to consider every orbit position and return a 7500xM rather than a single
%       percent success for each orbit position. Function "MixOrbitPositions" is used to
%       evaluate the output for multiple systems
%V4 - updated 2/8/2021
%   added unit converstion to km/s for DV1,DV2
%   updated burn time fraction in check transfer to 0.5
%   changed DV adjustment to use polynomial, changed DVadjust input to p1_adjust and
%       p2_adjust
%   changed maxflyby from a constant to polynomial, added input p_flyby






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





f_1 = exp(DV1/(g*-Isp1));
f_2 = exp(DV2/(g*-Isp2));

m2 = m3/f_2;
m_dump2 = m3 + m_prop2 - m2;


dt_burn = g*Isp1*m_burned1/T1;

        
DV1_ii = -Isp1*g*log(f_1) -Isp1*g*log(f_1)



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


