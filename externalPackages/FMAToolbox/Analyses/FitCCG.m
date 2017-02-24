function [dt1,dt2,Travel,index] = FitCCG(t,ccg,varargin)

%FitCCG - Fit dampened sinewave to the cross-correlogram of a pair of theta-modulated cells.
%
%  USAGE
%
%    [dt1,dt2,index] = FitCCG(t,ccg,<options>)
%
%    t              time bins
%    ccg            cross-correlogram
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'show'        display and plot results (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    dt1            non-parametric estimate of first theta-modulated peak
%    dt2            parametric estimate of first theta-modulated peak
%    index          theta modulation index (see Royer et al. 2010)

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
show = 'off';

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FitCCG">FitCCG</a>'' for details).');
end

if ~isdvector(t) || ~isdvector(ccg),
  error('Variables must be vectors (type ''help <a href="matlab:help FitCCG">FitCCG</a>'' for details).');
end
if length(t) ~= length(ccg),
  error('Variables are of different lengths (type ''help <a href="matlab:help FitCCG">FitCCG</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FitCCG">FitCCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'off','on'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help FitCCG">FitCCG</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help FitCCG">FitCCG</a>'' for details).']);
	end
end

global FitCCG_t FitCCG_ccg;
FitCCG_t = t;
FitCCG_ccg = ccg;

ccg = ccg(:);

% First estimate of dt: smooth the CCG and find the first peak for t>0
smoothed = Smooth(ccg,3);
maxima = IsExtremum([t smoothed]);
peaksPos = t(t>0&maxima);
peaksNeg = t(t<0&maxima);

if peaksPos(1) < abs(peaksNeg(size(peaksNeg,1))),
    dt1 = peaksPos(1);
else 
    dt1 = peaksNeg(size(peaksNeg,1));
end

% Try to find decent initial values
tau = 0.5;
mu = 0;
omega = 7; % near theta
b = mean(ccg);
phi = 2*pi*omega*dt1;
a = b;

beta = [a tau mu omega phi b];

% Compute best fit
[beta,~,code] = fminsearch(@Fit,beta);
a = abs(beta(1));
tau = -abs(beta(2));
mu = beta(3);
omega = beta(4);
phi = beta(5);
b = beta(6);

% Get the fit corresponding to the estimated parameters
[~,fit] = Fit(beta);

% Second estimate of dt: find the first peak of the model for t>0 or t<0
% (the closest to zero)
maxima = IsExtremum([t fit]);
peaksPos = t(t>0&maxima);
peaksNeg = t(t<0&maxima);

if peaksPos(1) < abs(peaksNeg(size(peaksNeg,1))),
    dt2 = peaksPos(1);
else 
    dt2 = peaksNeg(size(peaksNeg,1));
end

% Theta index
index = a/b;

% Find T (Geisler 2011) : Travel time
carrier = smooth(fit,5);
[~,I] = max(carrier);
Travel = t(I);

% Plot and display results
if strcmp(show,'on'),
	printf('a       = %.3f',a);
	printf('tau     = %.3f',tau);
	printf('mu      = %.3f',mu);
	printf('omega   = %.3f',omega);
	printf('phi     = %.3f',phi);
	printf('b       = %.3f',b);

	figure;
	hold on;
	bar(t,ccg);
	plot(t,fit,'r','linewidth',5);
end

% ------------------------------------------------------------------------------------------------
% Compute estimation error using a dampened sinewave model

function [error,model] = Fit(beta)

global FitCCG_t FitCCG_ccg;
t = FitCCG_t;
ccg = FitCCG_ccg;

a = beta(1);
tau = beta(2);
mu = beta(3);
omega = beta(4);
phi = beta(5);
b = beta(6);

% Fit a dampened sinewave to the data
model = (abs(a)*(sin(2*pi*omega*t + phi)+1) + b) .* exp(-abs(tau)*abs(t-mu));
% Compute summed square error
error = sum((ccg-model).^2);







%  options = optimset('Display','iter','MaxFunEvals',50,'MaxIter',50,'TolX',1e-6,'TolFun',1e-6,'PlotFcns',@optimplotfval);%'FunValCheck');
%  [beta,~,code] = fminsearch(@Fit,beta);%,options);
