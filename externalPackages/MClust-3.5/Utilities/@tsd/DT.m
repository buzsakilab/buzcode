function R = DT(tsa, tsflag)
%
% dt = tsd/DT(tsd, tsflag)
%	returns DT from tsa
% INPUTS
%      tsd 
%      tsflag: if 'ts' returns time in timestamps (default),
%              if 'sec' returns time in sec
%              if 'ms' returns time in ms
%
% ADR 
% version L2.0
% v 2.0: JCJ 3/3/2003 includes support for time units and robustness to
%        outliers.

if nargin == 2
    R = tsaTimes(tsa,tsflag,'d');
else
    R = nanmedian(diff(tsa.t));
end


