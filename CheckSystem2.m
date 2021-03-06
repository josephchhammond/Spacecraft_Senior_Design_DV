function [Success] = CheckSystem2(OrbitName,PosNum,Vsys,DV1sys,DV2sys,preposition_DV,tburn,mass_sys,R2sys,p1_adjust,p2_adjust,p_flyby,MainData,prop_scheme)
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
    DV1sys = DV1sys; %convert to km/s
    DV2sys = DV2sys; %convert to km/s
    Success = []; %create empty variable
    
    
    

    
    
    
    F1 = prop_scheme(1,1); %Departure: F=thrust [N], Isp=specific impulse [s]
    Isp1 = prop_scheme(1,3); %Departure: F=thrust [N], Isp=specific impulse [s]
    Isp2 = prop_scheme(2,3); %Arrival: Isp=specific impulse [s]
    
    Isp_dump1 = Isp1; %s ASSUMPTION
    Isp_dump2 = 300; %s ASSUMPTION

    m_stage = mass_sys(10);
    m_inert = mass_sys(8) - mass_sys(2) - mass_sys(3) - m_stage;
    m_prop1 = mass_sys(2);
    m_prop2 = mass_sys(3) - mass_sys(9);
    dtburn_dt_max = 0.5; % *ASSUMPTION ALERT*


    
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
        [pass_fail] = checktransfer(DV1data,DV2data,dtdata,R2data,Vsys,DV1sys,DV2sys,preposition_DV,tburn,R2sys,p1_adjust,p2_adjust,p_flyby,F1,Isp1,Isp2,Isp_dump1,Isp_dump2,m_inert,m_prop1,m_prop2,dtburn_dt_max,m_stage); %checks the ISO data agianst system data
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

function [pass_fail] = checktransfer(DV1data,DV2data,dtdata,R2data,Vsys,DV1sys,DV2sys,preposition_DV,tburn,R2sys,p1_adjust,p2_adjust,p_flyby,F1,Isp1,Isp2,Isp_dump1,Isp_dump2,m_inert,m_prop1,m_prop2,dtburn_dt_max,m_stage)
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





    
     g = 9.81; %m/s


% Convert dtdata into a 30x30 by extending it down each column
    dtdata_ = repmat(dtdata,30,1);

   
 
% DV2 (fourth burn) 
    % Reduce arrival requirements based on flyby velocity
    maxflyby = polyval(p_flyby,R2data);
    DV2adj_ = max(DV2data - maxflyby,0.001);
    
    % Find how many times we can cover the arrival requirements
    DV2_factor = DV2sys./DV2adj_/1000;
    
    % Find how much burntime we can allow for this system and actual DV
    % requirements for our arrival stage
        % Bad systems cant arrive under any circumstances and have unreal dtburn_dt < 0
        % Good systems can arrive under some amount of non-instantaenty and should have dtburn_dt 0-.5
        % Great systems are have extra fuel they dump to meet only the max burn time of dtburn_dt .5 
        

    % Set max burn time at the lower value of:
    %   a) Defined amount of 50% of transfer
    %   b) Max burn time from a departure system
    dt_burn_max = m_prop1*g*Isp1/F1;
    dtburn_dt_max_capable = dt_burn_max./dtdata_;
    check = dtburn_dt_max_capable < dtburn_dt_max;
    dtburn_dt_max = dtburn_dt_max + (dtburn_dt_max_capable - dtburn_dt_max) .* check;
    
    [dtburn_dt,DV2_factor_adj] = DV_adjustment2(DV2_factor,p2_adjust,dtburn_dt_max); 
    DV2_req = DV2_factor_adj.*DV2adj_; %Arrival DV specific to this transfer
    dtburn = dtdata_.*dtburn_dt;
    % Find fuel dumped
    f_2 = exp(DV2_req./(g*-Isp2));
    m_initial2 =m_inert./(f_2);
    m_burned2 = m_initial2 - m_inert;
    m_dump2 = m_prop2 - m_burned2;

% DV1 (third burn)
    % Based on allowable burn time, find burned fuel on departure stage and
    % how much is saved
    m_burned1 = dtburn.*F1./(g*Isp1);
    m_dump1 = m_prop1 - m_burned1;
    % Find how much departure DV we get, adjusted for non-instantaety
    f_1 = (m_inert + m_stage+ m_burned2)./(m_inert + m_stage + m_burned1 + m_burned2); %%%
    DV1_noninst = -Isp1.*g.*log(f_1);
    DV1_factor = DV_adjustment(1,p1_adjust,dtburn_dt,1);
    DV1_noninst_adj = DV1_noninst ./ DV1_factor;
    
% DV2 dump (second burn)
    % Find how much DV we get out of dumping extra arrival fuel before
    % departure
    m_post_dump2 = m_inert + m_burned1 + m_burned2 + m_stage;
    m_pre_dump2 = m_post_dump2 + m_dump2;
    f_dump2 = m_post_dump2./m_pre_dump2;
    DV_dump2 = -Isp_dump2*g.*log(f_dump2);
    
% DV1 dump (first burn)
    % Find how much DV we get out of dumping extra departure fuel before
    % departure
    m_post_dump1 = m_pre_dump2;
    m_pre_dump1 = m_post_dump1 + m_dump1;
    f_dump1 = m_post_dump1./m_pre_dump1;
    DV_dump1 = -Isp_dump1.*g.*log(f_dump1);
    
% Tally total departure DV output
    DV1_total = (DV_dump1 + DV_dump2 + DV1_noninst_adj - preposition_DV)/1000; %km/s
    
    
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



