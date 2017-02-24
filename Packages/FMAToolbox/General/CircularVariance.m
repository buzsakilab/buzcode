function [v,s] = CircularVariance(angles,dim)

%CircularVariance - Estimate circular variance and standard deviation.
%
%  USAGE
%
%    [v,s] = CircularVariance(angles,dim)
%
%    angles         angles in radians
%    dim            optional dimension along which the mean should be computed
%
%  OUTPUT
%
%    v              variance
%    s              standard deviation
%
%  SEE
%
%    See also CircularMean, CircularConfidenceIntervals, Concentration,
%    ConcentrationTest.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 2,
	dim = 1;
end

isradians(angles);

angles = exp(i*angles);
r_bar = abs(mean(angles,dim));
v = 1 - r_bar;
s = sqrt(-2*log(r_bar));
