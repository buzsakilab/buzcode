function [h,p,f] = FisherTest(samples1,samples2,alpha)

%FisherTest - Test if two groups of samples have equal variances.
%
% Fisher's test is used to test if two groups of samples have equal variances.
%
%  USAGE
%
%    [h,p,f] = FisherTest(samples1,samples2,alpha)
%
%    samples1       first group
%    samples2       second group
%    alpha          optional significance level (default = 0.05)
%
%  OUTPUT
%
%    h    test result (1 = reject null hypothesis, 0 = accept)
%    p    p value
%    f    test statistics

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if isempty(samples1) | isempty(samples2),
	h = 0;
	p = 0;
	f = 0;
	return
end

% Significance level
if nargin < 3,
	alpha = 0.05;
end

% Compute test
n1 = prod(size(samples1));
x1 = reshape(samples1,n1,1);
n2 = prod(size(samples2));
x2 = reshape(samples2,n2,1);
var1 = var(x1);
var2 = var(x2);
f = var1/var2;
df1 = n1-1;
df2 = n2-1;
if f > 1,
   p = 2*betainc(df2/(df2+df1*f),df2/2,df1/2);
else
   f = 1/f ;
   p = 2*betainc(df1/(df1+df2*f),df1/2,df2/2);
end
if p > 1,
   p = 2-p;
end

h = p < alpha;

if h,
	message = '+++ Variances are significantly different';
else
	message = '--- Variances are not significantly different';
end
message = ['F test: ' message ' (p='  num2str(p) ', F=' num2str(f) ', df1=' int2str(df1) ...
	', df2=' int2str(df2) ')'];
disp(message);
