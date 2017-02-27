function yred=rednoise(N,g,a);
% REDNOISE: A fast rednoise generator using filter.
%
% example: rednoise(N,g);
%
% Description: generates a rednoise series
% with zero process mean.
% note: statistical mean&variance will be different.
% 
% Inputs: n - desired length of time series
%         g - lag-1 autocorrelation
%         a - noise innovation variance parameter (optional, default=1)
%
% Example:
%   plot(rednoise(1000,.95))
%
% Aslak Grinsted 2006-2007

if nargin<3
    a=1;
end

if g==0
    yred=randn(N,1)*a;
    return
end
tau=ceil(-2/log(abs(g))); %2 x de-correlation time

yred=filter([1 0],[1;-g],randn(tau+N,1)*a);
yred=yred(tau+1:end);
