function [AverageCover,Spacing] = BestCombo(CoverageMatrix)
%BESTCOMBO Finds the best coverage percent and spacing for a given system with 2
%spacecraft. Will allow for spacing to be whatever it needs to be to return the maximum
%result, does not assume 180 degree spacing. If it is found this spacing is always best 
%or this function is too slow, then Joey's previous method should be used.
%
%Inputs:
%CoverageMatrix - nxn matrix, direct output of "MixOrbitPositions"
%
%Outputs:
%AverageCover - scalar value, average coverage percentage of the two spacecraft system
%Spacing      - sclalr value, number of positions on orbit from one spacecraft to the
%                   other
%
%Written by Colton Hook on 2/15/21

[n,~] = size(CoverageMatrix);
AverageCover_ = zeros(n,1);
for ii = 1:n
    running_sum = 0; %set running total to zero
    space = ii-1; %set spacing, start at zero
    for jj = 1:n
        spot1 = jj; %positon of S/C 1
        spot2 = mod(spot1+space-1,n)+1; %location of S/C 2, -1 and +1 are needed to keep spot2 on interval [1,n]
        running_sum = running_sum + CoverageMatrix(spot1,spot2);
    end
    AverageCover_(ii) = running_sum/n;
end
[AverageCover,Spacing] = max(AverageCover_);
Spacing = Spacing - 1;
end

