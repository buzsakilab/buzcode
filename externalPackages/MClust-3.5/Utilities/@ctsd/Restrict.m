function X = Restrict(X, t0, t1)

% ctsd/Restrict
% 	R = Restrict(tsa, t0, t1)
% 	R = Restrict(tsa, t0, t1)
% 	Returns a new tsa (ctsd) R so that D.Data is between 
%		timestamps t0 and t1, where t0 and t1 are in units
%
%   Assumes t has same units as D
%
% ADR 
% version L6.0
% v4.1 29 oct 1998 now can handle nargin=2
% v5.0 30 oct 1998 time dimension is always 1st dimension
% v5.1 19 jan 1998 now can handle t0 and t1 as arrays
% v6.0 JCJ 2/27/2003 includes support for time units
% v7.0 ADR 2007, sped up, reduced functionality to key factors

if nargin ~= 3
    error('Miscalled Restrict.');
end

if length(t0) ~= length(t1)
	error('t0 and t1 must be same length')
end

T = Range(X);
D = Data(X);

keep = false(size(T));
for iM = 1:length(t0)
	keep(T >= t0(iM) & T <= t1(iM)) = true;
end

T = T(keep);

% select along first dimension
shape = size(D);
D = reshape(D, [shape(1), prod(shape(2:end))]);
D = D(keep,:);
D = reshape(D, [size(D,1), prod(shape(2:end))]);

X = tsd(T, D);

