function [p1,p2] = NonInstantaneousLambert(Orbit)
% This function returns polynomials to be used to scale lambert solution
% maneuver DV requirements based on the orbit and "non-instantaenty"
% 
% PLEASE NOTE: Hardcoded in temporarily, function under maintenence!
%
% Inputs:
% Orbit = [perigee, apogee] in AU
% 
% Outputs:
% p1 = polynomials driving relationdship of dt_burn/dt to
% DV1_impulsive/DV1_nonimpulsive
% p2 = polynomials driving relationdship of dt_burn/dt to
% DV2_impulsive/DV2_nonimpulsive
% Written by: Joey Hammond

num = 1000;
t_linspaced = linspace(0,1,num);
for ii = 1:num
    dv1(ii) = .85*sqrt(t_linspaced(ii))+1 - t_linspaced(ii)/4;
    dv2(ii) = 1.05+.15*2*t_linspaced(ii);
end
p1 = polyfit(t_linspaced,dv1.^2,3);
p2 = polyfit(t_linspaced,dv2.^2,3);
end

