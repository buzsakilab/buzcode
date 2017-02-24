function R = EndTime(tsa, tsflag)
%
% R = tsd/EndTime(tsd, tsflag)
%	returns last timestamp covered by tsa
%
%      tsflag: if 'ts' returns time in timestamps (default),
%              if 'sec' returns time in sec
%              if 'ms' returns time in ms
%
%
%
% ADR 1998
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
    R = max(tsa.t);
end



