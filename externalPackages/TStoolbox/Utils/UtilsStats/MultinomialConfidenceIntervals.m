%MultinomialConfidenceIntervals - Simultaneous multinomial confidence intervals.
%
% Confidence intervals computed simultaneously on all cells using the method of
% Fitzpatrick and Scott (1987). This provides better estimates of the confidence
% intervals than the common (but erroneous) method where each cell is treated
% as an independent binomial variable.
%
%  USAGE
%
%    [p,boudaries] = MultinomialConfidenceIntervals(samples,alpha)
%
%    samples        list of numbers of samples in each cell
%    alpha          optional significance level (default = 0.05)
%
%    p              list of estimated cell probabilities
%    boundaries     confidence interval boudaries

% Copyright (C) 2004-2006 by Micha�l Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function [p,boundaries] = MultinomialConfidenceIntervals(samples,alpha)

% Parameter checking and default values

if nargin < 1,
  error('Incorrect number of parameters (type ''help MultinomialConfidenceIntervals'' for details).');
end

if nargin < 2
	alpha = 0.05;
end

% Process

n = sum(samples);
p = samples/n;

% There is a problem in the paper, so we cannot compute the intervals in the general case
% Thus, the following code is commented out

%  if alpha < 0.016,
%  	alpha = alpha/2;
%  elseif alpha < 0.15,
%  log(sqrt(2*pi)*(1-alpha/6))
%  	x = sqrt(-4/9*log(sqrt(2*pi)*(1-alpha/6)));
%  	alpha = 2*normcdf(x)
%  %  	alpha = 6*normcdf(3*norminv(alpha/2)/sqrt(8))-5;
%  else
%  	error('Cannot compute confidence intervals for this alpha level');
%  end
%  b = abs(norminv(alpha/2)/(2*sqrt(n)));

if alpha == 0.1,
	k = 1;
elseif alpha == 0.05,
	k = 1.13;
elseif alpha == 0.01,
	k = 1.40;
else
	error('Cannot compute confidence intervals for this alpha level');
end

b = k/sqrt(n);
boundaries = [max(p-b,0);min(p+b,1)];