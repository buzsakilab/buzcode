function kappa = Concentration(angles)

%Concentration - Estimate the concentration parameter for circular data.
%
% Estimate the concentration parameter for circular data assuming a Von Mises distribution.
% Uses the approximation described in "Statistical Analysis of Circular Data" (Fisher, p. 88).
%
%  USAGE
%
%    kappa = Concentration(angles)
%
%    angles         angles in radians
%
%  SEE
%
%    See also ConcentrationTest, CircularMean, CircularVariance, CircularConfidenceIntervals.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

isradians(angles);

n = length(angles);
angles = exp(i*angles);
r_bar = abs(mean(angles));

if r_bar < 0.53,
	kappa = 2*r_bar+r_bar^3+5*r_bar^5/6;
elseif r_bar < 0.85,
	kappa = -0.4+1.39*r_bar+0.43/(1-r_bar);
else
	kappa = 1/(r_bar^3-4*r_bar^2+3*r_bar);
end

% Correction for small samples
if n <= 15,
	if kappa < 2,
		kappa = max([kappa-2/(n*kappa) 0]);
	else
		kappa = (n-1)^3*kappa/(n^3+n);
	end
end
