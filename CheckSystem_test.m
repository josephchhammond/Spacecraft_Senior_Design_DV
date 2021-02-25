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
%Written by Colton Hook on 2/8/2021 - This is V4


%convert dtdata into a 30x30 by extending it down each column
    dtdata_ = repmat(dtdata,30,1);

%crit means criteria

    crit1 = R2data<R2sys; %we are capable than going further than this ISO
    
%adjust DV values
    maxflyby = polyval(p_flyby,R2data);
    DV2adj_ = max(DV2data - maxflyby,0);
    
    

    %While loop alculates fuel dumping and burn time adjustment on an orbit
    %by orbit basis
    error = 10^10;
    tburn1_ii = tburn1;
    g = 9.81; %m/s^2

    while error > 3600*24 %get on average, within a day of burntime accuracy
        %Adjust flyby mission REQUIRMENTS for given burntime
        DV1_req_ii = DV_adjustment(DV1data,p1_adjust,tburn1_ii,dtdata);
        DV2_req_ii = DV_adjustment(DV2adj_,p2_adjust,tburn1_ii,dtdata);
        %Dump arrival fuel based on mission requirements to find SYSTEM DV1 for a
        %given DV2 REQUIRMENT
        DV1_sys_ii = DV1_MAXsys + (DV1sys-DV1_MAXsys)/DV2sys * DV2_req_ii; %assuming a linear relationship between min and max fuel dump conditions, find required DV1 for this burn req

        % Find arrival prop mass dumped, if possible
        dm2_burned =  max(m2_burned - m2 * (1-exp(-DV2_req_ii*1000/(Isp2*g))),0);
        % Find departure prop mass dumped, if possible
        %dm1_burned = max(m1_burned - tburn1_ii*F1 / (g*Isp1),0); %kg
        dm1_burned =  max(m1_burned - m1 * (1-exp(-DV1_req_ii*1000/(Isp1*g))),0);
        % Find new starting mass FOR A GIVEN SYSTEM DV1 and DV2
        m1_adj = m1 - dm1_burned - dm2_burned;
        
        %Find new departure burn time FOR THE REQUIREMENT
        tburn1_adj = m1_adj*Isp1*g/F1 * (1-exp(-DV1_req_ii*1000/(Isp1*g))); %s
        %Find error and continue iteration
        error = mean(mean((tburn1_adj - tburn1_ii)));
        tburn1_ii = tburn1_adj;
    end

    fraction = 0.5; %change to higher fidelity model later *ASSUMPTION ALERT*
    crit2 = fraction*dtdata_ >= tburn1; %our burn is less than some fraction of the transfer
    crit3 = fraction*dtdata_ >= tburn2;

    
    DV1sys  = DV1_sys_ii;
    DV2sys = DV2sys;
    DV1req = DV1_req_ii;
    DV2req = DV2_req_ii;
    
    crit4 = DV2req <= DV2sys; %check system has enough fuel for arrival burn
    crit5 = DV1req <= DV1sys; %that a best case fuel dump for both stages is possible for departurne burn
    
%Check Volume
    crit6 = Vsys <= 1;
    
    pass_fail = crit1 & crit2 & crit3 & crit4 & crit5 & crit6; %which transfers met all 5 criteria
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

