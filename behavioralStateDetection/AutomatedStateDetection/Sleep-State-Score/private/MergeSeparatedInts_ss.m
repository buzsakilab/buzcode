function [ ints ] = MergeSeparatedInts( ints,minseparation )
%[ ints ] = MergeSeparatedInts( ints,minseparation ) merges close or
%overlapping intervals with a minimum separation.
%
%INPUT
%   ints            (n_ints x 2) matrix of interval start and end times
%   minseparation   merge ints separated by less than or equal to min
%
%DLevenstein 2016

%% Test
% ints = [0 3;2 4; 6 10; 10.5 13;  20 30; 19 21];
% minseparation = 0.5;


%%
%Sort Ints by their end times so all intervals end before the
%subsequent interval ends
[~,sort_ints] = sort(ints(:,2));
ints = ints(sort_ints,:);

%Find ints that start before a minimum separation after the end of the 
%previous interval
intseparation = ints(2:end,1)-ints(1:end-1,2);
smallsep = find(intseparation<=minseparation);

for s = flipud(smallsep)'   %Go backwards from last ints to first ones
    ints(s,2) = ints(s+1,2);
    ints(s,1) = min(ints(s+1,1),ints(s,1)); %Take the earlier start time between the ints
    ints(s+1,:) = [];
end

end

