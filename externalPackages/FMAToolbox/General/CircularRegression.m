function [beta,R2,beta_ts,R2_ts] = CircularRegression(x,angles,slope)

%CircularRegression - Non-parametric linear-circular regression.
%
%  Fit 'barber-pole' lines to the data. The model consists of three parallel lines
%  y = a.x + b + k*2pi (k in {-1,0,1}). The distance of a point to the model is
%  defined as the minimal distance to either of the three lines. The best-fit model
%  is computed using a least squared error approach.
%
%  Subsequently, reorganize angles around regression line (move by multiples of 2.pi)
%  and compute the more robust Theil–Sen estimators.
%
%  USAGE
%
%    [beta,R2,beta_ts,R2_ts] = CircularRegression(x,angles,slope)
%
%    x              linear independent variable
%    angles         dependent variable (in radians)
%    slope          optional search start value for slope (default = 0)
%
%  OUTPUT
%
%    beta           regression slope and intercept
%    R2             coefficient of determination
%    beta_ts        Theil–Sen slope and intercept
%    R2_ts          Theil–Sen coefficient of determination
%
%  SEE
%
%    See also CircularVariance, CircularConfidenceIntervals, Concentration,
%    ConcentrationTest.

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
end
if nargin < 3,
	slope = 0;
end
if ~isdvector(angles) || ~isdvector(x),
  error('Variables must be vectors (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
end
if length(angles) ~= length(x),
  error('Variables are of different lengths (type ''help <a href="matlab:help CircularRegression">CircularRegression</a>'' for details).');
end
isradians(angles);

x = x(:);
phi = angles(:);

beta = [];
R2 = [];
if isempty(angles) || isempty(x), return; end

% Store data in global variables (necessary for minimization)
global CircularRegression_x CircularRegression_phi;
CircularRegression_x = x;
CircularRegression_phi = phi;

% Minimize residual error (least squared error)
[beta,RSE] = fminsearch(@ResidualSquaredError,[slope 0]);

% Compute coefficient of determination
TSE = norm(phi-CircularMean(phi))^2;
R2 = 1-RSE/TSE;

% ------------------------------------------------------------------------------------------------
% Compute Theil–Sen estimator

if nargout <= 2, return; end
beta_ts = [nan nan];
R2_ts = nan;

% Move phi values by multiples of 2.pi to get them as close as possible to the regression line
d = phi-beta(1)*x+beta(2);
rd = sign(d).*floor(abs(d/(2*pi)));
phi = phi+rd*2*pi;

% Keep at most 500 random samples
n = length(x);
if n == 1, return; end
p = randperm(n);
n = min([n 500]);
x0 = x(p(1:n));
phi0 = phi(p(1:n));

% Compute pairwise slopes
slopes = nan(nchoosek(n,2),1);
k = 1;
for i = 1:n-1,
	for j = (i+1):n,
		slopes(k) = (phi0(i)-phi0(j))/(x0(i)-x0(j));
		k = k + 1;
	end
end

% Take median
beta_ts(1) = median(slopes);
beta_ts(2) = median(phi-beta_ts(1)*x);
% Compute coefficient of determination
RSE_ts = sum((phi-(beta_ts(1)*x+beta_ts(2))).^2);
R2_ts = 1-RSE_ts/TSE;

%  figure;hold on;plot([x;x],[phi;phi+2*pi],'.');
%  xm = mean([min(x);max(x)]);
%  ym = beta(1)*xm+beta(2);
%  ym_ts = beta_ts(1)*xm+beta_ts(2);
%  xLim = xlim;
%  PlotSlope(xm,ym,beta(1),diff(xLim),'k');
%  PlotSlope(xm,ym+2*pi,beta(1),diff(xLim),'k');
%  PlotSlope(xm,ym_ts,beta_ts(1),diff(xLim),'r');
%  PlotSlope(xm,ym_ts+2*pi,beta_ts(1),diff(xLim),'r');
%  xlim(xLim);


% --------------------------------------------------------------------------------------------------------------

function RSE = ResidualSquaredError(beta)

% Retrieve data
global CircularRegression_x CircularRegression_phi;
x = CircularRegression_x;
phi = CircularRegression_phi;

% Model parameters: slope and intercept
a = beta(1);
b = beta(2);

% Three lines
y1 = a*x+b;
y2 = y1+2*pi;
y3 = y1-2*pi;

% Squared error
d = min([(phi-y1).^2 (phi-y2).^2 (phi-y3).^2],[],2);
RSE = sum(d);
