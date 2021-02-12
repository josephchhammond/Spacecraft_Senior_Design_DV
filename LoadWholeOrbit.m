function [maindata] = LoadWholeOrbit(orbitname,numpos)
%LOADWHOLEORBIT Summary of this function goes here
%   Detailed explanation goes here

addpath(orbitname);
maindata = cell(1,numpos);
for i = 1:numpos
    filestring = "Results_" + i +".mat";
    maindata{1,i} = load(filestring);
end
rmpath(orbitname);
end

