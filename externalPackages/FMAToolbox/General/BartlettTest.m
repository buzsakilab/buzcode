function [h,p,T] = BartlettTest(data,alpha)

%BartlettTest - Test if k groups of samples have equal variances (homogeneity of variances).
%
%  This test assumes that all populations follow a Gaussian distribution.
%
%  USAGE
%
%    [h,p,t] = BartlettTest(data,alpha)
%
%    data           Nx2 matrix of (observation,group) pairs
%    alpha          optional significance level (default = 0.05)
%
%    h              test result (1 = reject null hypothesis, 0 = accept)
%    p              p value
%    t              test statistics

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Significance level
if nargin < 2,
	alpha = 0.05;
end

% Number of groups
k = max(data(:,2));
for i = 1:k,
	group = data(:,2) == i;
	% Number of samples in this group
	n(i) = sum(group);
	% Variance for this group
	s2(i) = var(data(group,1));
end
% Total number of samples
N = sum(n);
% Pooled variance
S2 = sum((n-1).*s2)/(N-k);

% Test statistics (biased)
t = (N-k)*log(S2)-sum((n-1).*log(s2));
% Bias correction
C = 1+(1/(3*(k-1)))*(sum(1./(n-1))-1/(N-k));
% Test statistics (unbiased)
T = t/C;

% Chi square at alpha with (k-1) degrees of freedom
p = 1 - chi2cdf(T,k-1);

h = p < alpha;

disp(['Bartlett test: Variances: ' num2str(s2)]);

if h,
	message = '+++ Two or more variances are significantly different';
else
	message = '--- Variances are not significantly different';
end
message = ['Bartlett test: ' message ' (p='  num2str(p) ', T=' num2str(T) ', N=' int2str(N) ')'];
disp(message);
