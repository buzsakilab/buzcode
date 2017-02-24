function R = EndTime(tsa, tsflag)
%
% T1 = ctsd/EndTime(tsd, tsflag)
%
% INPUTS
%      tsd 
%      tsflag: if 'ts' returns time in timestamps (default),
%              if 'sec' returns time in sec
%              if 'ms' returns time in ms
% ADR 1998
%
% version L5.0
% v4.1 Fixed bug: if width of data was longer than time gave back wrong answer
% v5.0: JCJ 12/19/2002 includes support for time units 
%
% Status: IN PROGRESS
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


if nargin == 2
    R = tsaTimes(tsa,tsflag,'e');
else
    R = tsa.t0 + tsa.dt * (size(tsa.data,1)-1);
end


