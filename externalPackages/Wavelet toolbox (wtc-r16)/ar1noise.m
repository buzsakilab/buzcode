  function X=ar1noise(n,c,g,a)
% AR1NOISE - Generate zero-mean red noise.
% Syntax: X=ar1noise(n,c,g,a);
%
% Inputs: n - desired length of time series
%         c - number of time series to generate
%         g - lag-1 autocorrelation
%         a - noise innovation variance parameter
% 
% Output: X - n by c matrix of red noise series.
%
% Written by Eric Breitenberger.      Version 1/21/96
% Please send comments and suggestions to eric@gi.alaska.edu       
%

% Comment out this line if you want the
% same 'random' realization each time:
randn('state',fix(sum(100*clock*(randn+4))));

X=zeros(n,c);
X(1,:)=sqrt(a^2/(1-g^2))*randn(1,c);
z=a*randn(n,c);

for i=2:n
  X(i,:)=g*X(i-1,:)+z(i,:);
end

% Center the surrogates:
X=X-ones(n,1)*mean(X);

