function R = DT(tsa,tsflag)
%
% dt = ctsd/DT(tsd, tsflag)
%	returns DT from tsa
% INPUTS
%      tsd 
%      tsflag: if 'ts' returns time in timestamps (default),
%              if 'sec' returns time in sec
%              if 'ms' returns time in ms
%
% ADR 
% version L5.0
% v 5.0: JCJ 12/19/2002 includes support for time units 
% 
% Status: IN PROGRESS
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

if nargin == 2
    R = tsaTimes(tsa,tsflag,'d');
else
    R = tsa.dt;
end


