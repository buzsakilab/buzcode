function [ INT ] = IDXtoINT_ss( IDX ,numstates)
%IDXtoINT(IDX) Converts state indices to state on/offsets
%
%INPUT
%   IDX:    [t x 1] vector of state indices, where states are identified by
%           integers starting from 1. Times with IDX 0 will not be counted
%           in any interval INT
%   numstates (optional)  number of interval types (for use
%
%OUTPUT
%   INT:    {nstates} cell array of intervals - start and end times
%
%DLevenstein 2015-16
%%

if ~exist('numstates','var')
    numstates = max(IDX);
end

states = 1:numstates;

if isrow(IDX)
    IDX = IDX';
end

IDX = [0; IDX; 0];
for ss = 1:numstates
    statetimes = IDX==states(ss);
    INT{ss} = [find(diff(statetimes)==1) find(diff(statetimes)==-1)-1];
end



end

