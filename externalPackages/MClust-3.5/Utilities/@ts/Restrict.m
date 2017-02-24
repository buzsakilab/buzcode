function X = Restrict(X, t0, t1)

% 	R = Restrict(X, t0, t1)
% 	Returns a new ts R containing only samples between 
%		timestamps t0 and t1, where t0 and t1 are in units
%   NOTE: assumes units provided are the same as the units in the tsd
%   NOTE: D will be returned in its original units!
%
% ADR 2000
% version 6.0
% v4.1 29 oct 1998 now can handle nargin=2
% v4.2 23 jan 1999 now can handle t0 and t1 as arrays
% v5.0 JCJ 2/27/2003 includes support for time units
% v6.0 ADR 5/Jan/2007 removed units support, sped up 

if nargin ~= 3
    error('Miscalled Restrict.');
end

if length(t0) ~= length(t1)
	error('t0 and t1 must be same length')
end

keep = false(size(X.t));
for iM = 1:length(t0)
	keep(X.t >= t0(iM) & X.t <= t1(iM)) = true;
end

X.t = X.t(keep);