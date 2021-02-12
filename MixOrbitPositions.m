function [PercentCoverage] = MixOrbitPositions(SuccessList1,SuccessList2)
%MIXORBITPOSITIONS Looks at two orbits and returns the percent coverage if one
%spacecraft is placed in combinations of each position on both orbits.
%
%Inputs:
%SuccessList1 - NxM matrix of if the system can reach each ISO where each of the N rows is
%               an ISO and the M columns are orbit locations
%SuccessList2 - NxP matrix of if the system can reach each ISO where each of the N rows is
%               an ISO and the P columns are orbit locations
%note: N must match between inputs. M and P can be different
%
%Outputs:
%PercentCoverage - matrix of size PxM percent of ISOs covered by the system at the two orbit locations specifed                
%
%Written by Colton Hook on 2/8/2021 - This is V2

M = length(SuccessList1(1,:));
P = length(SuccessList2(1,:));

PercentCoverage = zeros(M,P); %initalize square matrix YxY, one location for each orbit
for ii = 1:M
    for jj = 1:P
        spot1 = SuccessList1(:,ii);
        spot2 = SuccessList2(:,jj);
        n = length(spot1); %how many ISOs are there
        combined = spot1|spot2; %each row is an ISO, is 1 if *at least 1* of the spacecraft can reach it, 0 otherwise
        PercentCoverage(ii,jj) = 100 * sum(combined)/n;
    end
end
end

%---------------------------Change Log---------------------------
%V1 created 2/5/2021
%V2 updated 2/8/2021
%   added input for second SuccessList, no checks two different success lists, use the
%       same one for both inputs if only checking one orbit
%   removed ability to check a range or positions, now checks all be default
