function [h,U2] = WatsonU2Test(group1,group2,alpha)

%WatsonU2Test - Test if two samples (circular data) have different means / variances.
%
%  This non-parametric test assumes the data comes from a continuous distribution.
%
%  USAGE
%
%    [h,U2] = WatsonU2Test(group1,group2,alpha)
%
%    group1         angles in radians for group 1
%    group2         angles in radians for group 2
%    alpha          optional significance level (default = 0.05)
%
%    h              test result (1 = reject null hypothesis, 0 = accept)
%    U2             test statistics
%

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

isradians(group1);
isradians(group2);

% Significance level
if nargin < 3,
	alpha = 0.05;
end

group1 = group1(:);
group2 = group2(:);

% Number of angles in each group
n1 = length(group1);
n2 = length(group2);
if n1 < 5 || n2 < 5,
	error('This test requires a minimum of 5 angles per group');
end
N = n1+n2;

% Build an Nx2 matrix:
%  column 1: all data in group 1 and group 2 in ascending order
%  column 2: 0 if data comes from group 1, 1 otherwise
data = [group1 zeros(size(group1));group2 ones(size(group2))];
sorted = sortrows(data);
% Build i and j, the respective cumulative functions for groups 1 and 2
i = zeros(N,1);
j = zeros(N,1);
ingroup1 = sorted(:,2)==0;
i(ingroup1) = 1;
i = cumsum(i)/n1;
j(~ingroup1) = 1;
j = cumsum(j)/n2;
% Compute dk
dk = i-j;
% Test statistic
U2 = n1*n2/N^2*(sum(dk.^2)-sum(dk)^2/N);

% Perform test

% Number of elements in each group, modified to fit table entries
k1 = min([n1 n2]);
k2 = max([n1 n2]);
if k1 > 10, k1 = inf; end
if k2 > 12, k2 = inf; end
% Critical values (from Kanji, 100 statistical tests, 1999)
m1 = [5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 8 9 9 9 9 10 10 10 inf];
m2 = [5 6 7 8 9 10 11 12 6 7 8 9 10 11 12 7 8 9 10 11 12 8 9 10 11 12 9 10 11 12 10 11 12 inf];
if alpha == 0.01,
	a = [nan nan nan nan 0.28 0.289 0.297 0.261 nan 0.282 0.298 0.262 0.248 0.262 0.259 0.304 0.272 0.255 0.262 0.253 0.252 0.250 0.258 0.249 0.252 0.252 0.266 0.254 0.255 0.254 0.255 0.255 0.255 0.268];
elseif alpha == 0.05,
	a = [0.225 0.242 0.2 0.215 0.191 0.196 0.19 0.186 0.206 0.194 0.196 0.193 0.19 0.187 0.183 0.199 0.182 0.182 0.187 0.184 0.186 0.184 0.186 0.185 0.184 0.185 0.187 0.186 0.185 0.185 0.185 0.186 0.185 0.187];
else
	error('This test is implemented for alpha levels of 0.05 and 0.01 only.');
end
% Find critical value
i1 = k1 == m1;
i2 = k2 == m2;
value = a(i1&i2);
if isnan(value),
	error(['The test cannot be computed for n1=' int2str(n1) ', n2=' int2str(n2) 'and alpha=' num2str(alpha) '.']);
end

% Test result
h = U2 >= value;

disp(['Watson U2 test: Means:       m1=' num2str(CircularMean(group1)) ', m2=' num2str(CircularMean(group2))]);
disp(['Watson U2 test: Variances:   v1=' num2str(CircularVariance(group1)) ', v2=' num2str(CircularVariance(group2))]);

if h,
	message = '+++ The two groups have significantly different means/variances';
else
	message = '--- The two groups do not have significantly different means/variances';
end
message = ['Watson U2 test: ' message ' (n1='  int2str(n1) ', n2=' int2str(n2) ', U2=' num2str(U2) ', critical=' num2str(value) ')'];
disp(message);
