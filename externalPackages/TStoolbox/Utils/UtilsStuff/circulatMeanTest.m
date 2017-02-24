function Y = circularMeanTest(t1,t2,varargin)

% USAGE
%
%    [H,P] = circularMeanTest(t1,t2)
%	performs a on parametric test on the mean where t1 and t2 are two populations of circular data.
%    H and P are the boolean value of success (at P<0.05 if not specified) and P is the p-value.
% 
% From Watson 1983, in Fischer p116
% adapted to matlab by Adrien Peyrache, 2007

n = zeros(2,1);
m = zeros(2,1);

n(1) = length(t1);
n(2) = length(t2);
N = sum(n);

[m(1), Rmean, d(1), pval] = CircularMean(t1);
[m(2), Rmean, d(2), pval] = CircularMean(t2);

C = n'*cos(m);
S = n'*sin(m);

R = sqrt(C^2+S^2);
d = n'*d/N;

Y = 2*(N-R)/d;


