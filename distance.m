function [dist, index] = mydistance(current, Target)
%This function calculates the distance between any two cartesian 
%coordinates.

num = size(Target, 1);

[dist, index] = min(sqrt(sum((repmat(current, num, 1) - Target) .^ 2, 2)));

