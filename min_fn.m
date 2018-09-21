function min_point_index = min_fn(OPEN, OPEN_COUNT, xTarget, yTarget)
%Function to return the Node with minimum fn
% This function takes the list OPEN as its input and returns the index of the
% node that has the least cost
%
%   Copyright 2009-2010 The MathWorks, Inc.

 OPEN_valid = OPEN(:, 1) == 1;
 
 % whether reach the target
 Target = [xTarget, yTarget];
 Istarget = ismember(OPEN(:, 2:3), Target, 'rows');
 Istarget = Istarget & OPEN_valid;
 if sum(Istarget)
     % one of the successors is the goal node so send this node
     min_point_index = find(Istarget);
 end
 
 %Send the index of the smallest node
 if sum(OPEN_valid)
     [min_fn, min_index] = min(OPEN(OPEN_valid, 7));
     valid_index = find(OPEN_valid);
     min_point_index = valid_index(min_index);
 else
     min_point_index = -1;
 end