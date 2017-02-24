function R = TimePoints(tsa, tsflag)

% Returns TSD timestamps, calls "Range"
% 
%  USAGE:
%  R = timePoints(tsa)
%  R = timePoints(tsa, tsflag)
%  
%  OPTIONS:
%  tsflag - if 'ts' returns time in timestamps (default),
%  	    if 'sec' returns time in sec
%  	    if 'ms' returns time in ms


% Brendon Watson 
% version 1.0
% November 27, 2013

if nargin == 2;
    R = Range(tsa, tsflag);  
elseif nargin == 1;
    R = Range(tsa);
end