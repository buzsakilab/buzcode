function R = StartTime(tsa, tsflag)

% t = tsd/StartTime(tsd, tsflag)
%	returns first timestamp covered by tsa
% returns last timestamp in TS
%      tsflag: if 'ts' returns time in timestamps 
%              if 'sec' returns time in sec (default)
%              if 'ms' returns time in ms
%
%
% ADR
% version L5.0
% v5.0: JCJ 12/19/2002 includes support for time units 
%
% Status: IN PROGRESS
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


if nargin == 2
    R = tsaTimes(tsa,tsflag,'e');
else
    R = min(tsa.t);
end



