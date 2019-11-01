function [out, segmentLength] = suprathresh(V, thresh)
% This subroutine determines suprathreshold values
%%% INPUTS %%%
% V:        a Nx1 (checked) vector of values
% thres:    a scalar threshold
%
%%% OUTPUTS %%%
% out:      a Nx2 array with the 1st column indicating first suprathreshold
%           values and the 2nd column the last for each block
% segmentLength:    a Nx1 indicating the number of samples within each
%           suprathreshold block

% check for column vector
if size(V, 1) < size(V, 2)
    V = V';
end

% evaluate threshold
V = V >= thresh;
V = [0; V; 0];
d = diff(V);

% set output
out           = [find(d == 1), find(d == -1) - 1];
segmentLength = diff(out, [], 2) + 1;
