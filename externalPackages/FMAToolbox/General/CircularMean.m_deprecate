function [m,r] = CircularMean(angles,dim)

%CircularMean - Estimate the circular mean.
%
%  USAGE
%
%    [m,r] = CircularMean(angles,dim)
%
%    angles         angles in radians
%    dim            optional dimension along which the mean should be computed
%
%  OUTPUT
%
%    m              mean angle
%    r              resultant length
%
%  SEE
%
%    See also CircularVariance, CircularConfidenceIntervals, Concentration,
%    ConcentrationTest.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if isempty(angles), m = []; return; end

if nargin < 2,
	dim = 1;
end

isradians(angles);

angles = exp(i*angles);
r = abs(nanmean(angles,dim));
m = angle(nanmean(angles,dim));
