function [X, xi] = runavg(x, n, EXCLUDE)
% [X, xi] = runavg(x, n)  [n should be odd]
%
% X is the running average of x using an n-element window around each value
% in x (except for (n-1)/2 elements at beginning and end). If EXCLUDE flag
% is true, the average is taken around each point, with that point excluded
% from the average. xi gives the indices in x corresponding to the values
% in X.
%
% Created by Erik Schomburg, January 2012

if nargin < 3
    EXCLUDE = 0;
elseif (EXCLUDE > 0) && ~mod(EXCLUDE,2)
    error('If exclude arg is nonzero, it should be odd');
end

if ~mod(n,2)
    error('n arg should be odd');
end

if EXCLUDE > 0
    X = filter([ones(1,(n-EXCLUDE)/2) zeros(1,EXCLUDE) ones(1,(n-EXCLUDE)/2)]/(n-EXCLUDE),1,x);
else
    X = filter(ones(1,n)/n,1,x);
end
X = X(n:end);
xi = (1+(n-1)/2):(length(x)-(n-1)/2);
