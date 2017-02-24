function X = Mask(X, posneg, t0, t1)
%
% mTSD = ctsd/Mask(tsd, posneg, t0, t1)
%
% INPUTS:
%    X = ctsd object
%    posneg = 1=NaN times inside trial pairs, 0=NaN times outside trial pairs
%    t0 = sets of start times
%    t1 = sets of end times
%     NOTE: must be SAME units as tsd!
%
% OUTPUTS:
%    mtsd = tsd object with times (not) in TrialPairs set to NaN
%
% ADR 1998
% version L6.0 
% ADR Jan/2007 - incorporated with posneg
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

% Unwrap trial pairs

if length(t0) ~= length(t1)
	error('t0 and t1 must be same length')
end

T = Range(X);
D = Data(X);

keep = false(size(T));
for iM = 1:length(t0)
	keep(T >= t0(iM) & T <= t1(iM)) = true;
end

% select along first dimension
shape = size(D);
D = reshape(D, [shape(1), prod(shape(2:end))]);
if posneg
	D(keep,:) = NaN;
else %not posneg
	D(~keep,:) = NaN;
end	
D = reshape(D, [size(D,1), prod(shape(2:end))]);

X = tsd(T, D);
