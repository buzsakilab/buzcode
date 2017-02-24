function R = Range(tsa, tsflag)

% tsd/Range
% 
%  R = Range(tsa)
%  R = Range(tsa, 'sec')
%  R = Range(tsa, 'ts')
%  R = Range(tsa. 'all_ts')
%
%  returns range covered by tsa
%      tsflag: if 'ts' returns time in timestamps
%              if 'sec' returns time in sec 
%              if 'sec0' returns time in sec counting from 0
%              if 'ms' returns time in ms
%
% ADR 2000
% version L6.0
% v4.1 28 oct 1998 flag no longer optional.
% v 5.0: set up to expect inputs with ts in SEC (as per new ADRLAB policy)
% v 6.0: JCJ 12/19/2002 includes support for time units 
%
% Status: IN PROGRESS 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


if nargin == 2
    R = tsaTimes(tsa,tsflag,'a');
else
    R = tsa.t;
end


