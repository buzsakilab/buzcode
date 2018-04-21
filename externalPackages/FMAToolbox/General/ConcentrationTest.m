function [h,p] = ConcentrationTest(angles,group,alpha,nRandomizations)

%ConcentrationTest - Test homogeneity of concentration parameters.
%
% Test whether the concentration parameters for N groups of circular data are equal.
% The data are assumed to have Von Mises distributions. If the median concentration
% is inferior to 1, a randomization procedure is used.
% See "Statistical Analysis of Circular Data" (Fisher, p. 131).
%
%  USAGE
%
%    [h,p] = ConcentrationTest(angles,group,alpha,nRandomizations)
%
%    angles         angles in radians
%    group          group number for each sample
%    alpha          optional significance level (default = 0.05)
%    n              optional number of randomizations (default = 1000)
%
%  SEE
%
%    See also Concentration, CircularMean, CircularVariance, CircularConfidenceIntervals.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global n N r;

if nargin < 3,
	alpha = 0.05;
end
if nargin < 4,
	nRandomizations = 1000;
end
isradians(angles);

% Number of groups (r), number of samples in each group (n), total number of samples (N)
r = max(group);
N = length(angles);
for i = 1:r,
	n(i) = sum(group==i);
end
if min(n) < 10,
	error('This test requires a minimum of 10 samples per group.');
end

% Median concentration parameter (kappa)
for i = 1:r,
	k(i) = Concentration(angles(group==i));
end
kappa = median(k);

% Circular mean for each group (mu)
for i = 1:r,
	mu(i) = CircularMean(angles(group==i));
end

if kappa >= 1,
	% Statistics for the data
	fr = ComputeStatistics(angles,group,mu);
	% p level
	p = 1 - fcdf(fr,r-1,N-r);
	p = 2 * min(p, 1-p);
else
	% Center all samples
	for i = 1:r,
		angles(group==i) = angles(group==i)-mu(i);
	end
	% Compute the statistics for the data
	fr = ComputeStatistics(angles,group,mu);
	% Compute the statistics for each randomization
	for i = 1:nRandomizations,
		x = randperm(N);
		group = group(x);
		for j = 1:r,
			mu(j) = CircularMean(angles(group==j));
		end
		fr_rand(i) = ComputeStatistics(angles,group,mu);
	end
	% Sort statistics
	fr_rand = sort(fr_rand);
	% Where does the statistics for the original data fall?
	m = find(fr_rand>=fr);
	% p level
	if isempty(m),
		p = 0;
	else
		m0 = sum(m==m(1));
		m = m(1);
		if m0 == 1,
			p = (nRandomizations-m+1)/nRandomizations;
		else
			p = (nRandomizations-m)/nRandomizations+m0/2;
		end
	end
end

h = p < alpha;

message = ['Concentration test: median kappa = ' num2str(kappa) ' ('];
for i = 1:r,
	message = [message num2str(k(i)) ' '];
end
message(end) = ')';
disp(message);

if h,
	message = '+++ Two or more concentration parameters are significantly different';
else
	message = '--- Concentration parameters are not significantly different';
end
message = ['Concentration test: ' message ' (p='  num2str(p) ', fr=' num2str(fr) ', N=' int2str(N) ')'];
disp(message);


function fr = ComputeStatistics(angles,group,mu)

global n N r;

for i = 1:r,
	a = abs(sin(angles(group==i)-mu(i)));
	d(i) = sum(a)/n(i);
	s(i) = sum((a-d(i)).^2);
end

D = sum(n.*d)/N;

fr = (N-r)*sum(n.*(d-D).^2)/((r-1)*sum(s));
