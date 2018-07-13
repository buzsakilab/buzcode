function [ window ] = bz_RandomWindowInIntervals( intervals,winsize )
%[ window ] = RandomWindowInInterval( interval,winsize ) returns a random
%time window from within an interval set.
%
%INPUTS
%   interval     [N x 2] start stop pairs
%   winsize      duration of the window
%
%OUTPUTS
%   window       random time window from within the interval
%
%DLevenstein 2018
%%
%Intervals of possible start times
if isequal(size(intervals),[2 1])
    intervals = intervals';
end

intervals(:,2) = intervals(:,2)-winsize;
intervals(diff(intervals,1,2)<=0,:)=[];

maxt = max(intervals(:));

possiblestarttimes = 1:1:maxt;
possiblestarttimes = Restrict(double(possiblestarttimes),double(intervals));

starttime = randsample(possiblestarttimes,1);

window = starttime+[0 winsize];


end

