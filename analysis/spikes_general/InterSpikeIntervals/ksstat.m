function ks = ksstat(x,y)
% ks = ksstat(x,y)
% Computes the two-sample Kolmogorov-Smirnov statistic.
% x & y are two independent groups.
% ks is the maximum of the absolute differences between the two empirical
% cumulative distribution functions (CDF).
%
% No p value is returned.
%
% Based on Rand Wilcox's ks & ecdf R functions
% http://dornsife.usc.edu/labs/rwilcox/software/
%
% The Matlab function kstest2 returns the same results.

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% remove NaNs & reformat
x=x(~isnan(x)); x=x(:); Nx = numel(x);
y=y(~isnan(y)); y=y(:); Ny = numel(y);

all = sort([x;y]); % pool and sort the observations
diff = zeros(size(all)); 
for e = 1:numel(all)
    diff(e) = abs( length(x(x<=all(e)))/Nx - length(y(y<=all(e)))/Ny );
end    
ks = max(diff);