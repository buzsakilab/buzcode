% Conduct a binary search
% 
% binSearch: Conducts a binary search and continues until the score is
% close enough to the target. An upper limit should be set on the number of
% iterations
%
%     [SearchPoint, Step, ReLoop] = iosr.auditory.binSearch(Score, ...
%       TargetScore, ...
%       CloseEnough, ...
%       NextSearchPointa, ...
%       Step, ...
%       LoopCount, ...
%       MaxLoops);
% end
%
% inputs:
% - Score: the test value
% - Target: the 'finished' test value
% - CloseEnough: the distance from the Target at which it is acceptable to
% discontinue the binary search
% - SearchPoint: the next work input value to try
% - Step: the distance which the SearchPoint can move by
% - LoopCount: a counter for the number of iterations completed so far
% - MaxLoops: The upper limit on the number of iterations allowed
%
% outputs:
% - SearchPoint: the next work input value you should try
% - Step: the distance which can be stepped on the NEXT iteration. You
%   should store this value.
% - ReLoop: If this flag is set to 0, the upper level while loop will be
%   terminated
%
%
% example: you wish to use an audibility model to find the level at which
% you can be 50% confident of detecting the signal (-+ 1%). You choose to
% start with an input level of 20dB and take no more than 10 steps of
% 40,20,10 dB etc.
%
% You would implement this in the following way:
%
% TargetPercent = 0.5;
% CloseEnough = 0.01;
% Signal Level = 40;
% Step = 40;
% MaxLoops = 10;
% LoopCount = 0;
%
% while ReLoop = 1
%
%   LoopCount = LoopCount + 1;
%
%   Percent = RunAudibilityModel(SignalLevel)
%
%   [SignalLevel Step ReLoop] = iosr.auditory.binSearch(Percent, ...
%       TargetPercent, ...
%       CloseEnough, ...
%       SignalLevel, ...
%       Step, ...
%       LoopCount, ...
%       MaxLoops);
%
% end
% 

%   Copyright 2016 University of Surrey.

function [SearchPoint, Step, ReLoop] = binSearch(Score, Target, CloseEnough, SearchPoint, Step, LoopCount, MaxLoops)

% input tests
% test input types
assert(isnumeric(Score) ... 
    &isnumeric(Target) ...
    &isnumeric(CloseEnough) ...
    &isnumeric(SearchPoint) ...
    &isnumeric(Step) ...
    &isnumeric(LoopCount) ...
    &isnumeric(MaxLoops),'input arguments must be numeric!');


% Define next search point
if Score > Target
    SearchPoint = SearchPoint - Step;
else
    SearchPoint = SearchPoint + Step;
end

% half the distance of the next step
Step = Step / 2;

% if we got close enough to our taget, or it was the last iteration
if (abs((Score-Target))<CloseEnough)||(LoopCount>=MaxLoops)
    ReLoop = 0;
else
    ReLoop = 1;
end

end
