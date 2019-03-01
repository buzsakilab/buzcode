function [ window ] = bz_RandomWindowInIntervals( intervals,winsize,nwins )
%[ window ] = bz_RandomWindowInIntervals( interval,winsize,nwins ) returns a random
%time window from within an interval set.
%
%INPUTS
%   interval     [N x 2] start stop pairs
%   winsize      duration of the window
%   nwins        number of random windows to get
%
%OUTPUTS
%   window       random time window(s) from within the interval
%
%DLevenstein 2018
%%
if exist('nwins','var') 
    for nn = 1:nwins
        window(nn,:) = bz_RandomWindowInIntervals(intervals,winsize);
        %Should exclude the window from intervals here...
        [~,startin] = InIntervals(window(nn,1),intervals);
        [~,endin] = InIntervals(window(nn,2),intervals);
        newint_pre = [intervals(startin,1) window(nn,1)];
        newint_post = [window(nn,2) intervals(endin,2)];
        intervals(startin,:) = [];
        intervals = [intervals; newint_pre; newint_post];
    end
    return
end

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

