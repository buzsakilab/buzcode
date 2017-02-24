function delay = PETHTransition(psth,timeBins,varargin)

%PETHTransition - Find a transition point in a peri-event time histogram.
%
% Find a transition point in a peri-stimulus time histogram using the maximum
% likelihood or mean square error algorithm of Friedman and Priebe (1998), or
% piecewise linear fit.
%
%  USAGE
%
%    delay = PETHTransition(peth,timeBins,<options>)
%
%    peth           PETH (computed using <a href="matlab:help SyncHist">SyncHist</a>)
%    timeBins       time bins
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'method'      either 'ml' (maximum likelihood, default), 'ls' (least
%                   squares) or 'pl' (piecewise linear)
%     'show'        either 'on' (default) or 'off'
%    =========================================================================
%
%  EXAMPLE
%
%    [raster,indices] = Sync(spikes,stimuli);     % compute spike raster data
%    figure;PlotSync(raster,indices);             % plot spike raster
%    [s,t] = SyncHist(raster,indices);            % compute PETH
%    delay = PETHTransition(s,t);                 % find transition point
%
%  SEE
%
%    See also Sync, SyncHist, SyncMap, PlotSync.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
method = 'ml';
show = false;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PETHTransition">PETHTransition</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help PETHTransition">PETHTransition</a>'' for details).']);
  end
  switch(lower(varargin{i})),

    case 'method',
		method = lower(varargin{i+1});
      if ~isstring_FMAT(method,'ml','ls','pl'),
        error('Incorrect value for property ''method'' (type ''help <a href="matlab:help PETHTransition">PETHTransition</a>'' for details).');
      end

    case {'show','plot'},
    	show = lower(varargin{i+1});
    	if ~isstring_FMAT(show,'on','off'),
        error('Incorrect value for property ''show'' (type ''help <a href="matlab:help PETHTransition">PETHTransition</a>'' for details).');
      end

    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PETHTransition">PETHTransition</a>'' for details).']);

  end
end

if strcmp(method,'ml'),

	% Maximum likelihood

	f = psth;
	n = length(f);
	for i = 1:n,
		F(i) = sum(log(1:max([1 f(i)])));
	end
	for theta = 1:n-1,
		t1 = 1:theta;t1 = t1';
		t2 = theta+1:n;t2 = t2';
		lambda1 = sum(f(t1))/(theta+1);
		lambda2 = sum(f(t2))/(n-theta);
		likelihood(theta) =  ...
			- lambda1 * (theta+1) ...
			+ log(lambda1) * sum(f(t1)) ...
			- sum(F(t1)) ...
			- lambda2 * (n-theta) ...
			+ log(lambda2) * sum(f(t2)) ...
			- sum(F(t2));
	end
	likelihood(end+1) = likelihood(end); % dummy value
	d = find(likelihood == max(likelihood));
	d = d(1);
	delay = timeBins(d);

	if show,
		fig = figure;
		set(fig,'name','PETH Transition - Maximum Likelihood','number','off');
		subplot(2,1,1);hold on;
		bar(timeBins,psth);
		plot([delay delay],ylim,'k','linestyle','--');
		title(['delay = ' num2str(delay) ' s']);
		ylabel('Occurrences');
		subplot(2,1,2);hold on;
		plot(timeBins,likelihood,'Color',[1 0 0],'Marker','none','LineStyle','-');
		plot([delay delay],ylim,'k','linestyle','--');
		xlabel('Time (s)');
		ylabel('Likelihood');
	end

elseif strcmp(method,'ls'),

	% Least square error

	f0 = psth;
	n0 = length(f0);
	t0 = 1:n0;t0 = t0';
	for i = 1:n0,
		F0(i) = sum(f0(1:i-1));
	end
	F0 = F0';
	start = 1;
	stop = length(f0);
	f = f0(start:stop);
	n = length(f);
	for i = 1:length(f),
		F(i) = sum(f(1:i-1));
	end
	F = F';
	t = 1:n;t = t';
	squaredSum(1:n) = NaN;
	for theta = 2:n-1,
		t1 = 1:theta; t1 = t1';
		t2 = theta+1:n; t2 = t2';
		lambda1(theta) = ...
			(...
			- sum(t1 .* F(t1)) ...
			- theta * sum(F(t2)) ...
			+ theta * sum(theta-t2) * sum(F(t2) .* (theta-t2)) / sum((theta-t2).^2) ...
			) / (...
			+ theta^2 * sum(theta-t2) * sum(theta-t2) / sum((theta-t2).^2) ...
			- sum(t1.^2) ...
			- theta^2 * (n-theta) ...
			);
		lambda2(theta) = ...
			(...
			+ lambda1(theta) * theta * sum(theta-t2) ...
			- sum(F(t2).*(theta-t2)) ...
			) ...
			/ sum((theta-t2).^2);
		squaredSum(theta) =  ...
			+ sum((F(t1)-lambda1(theta)*t1).^2) ...
			+ sum((F(t2)-lambda1(theta)*theta+lambda2(theta)*(theta-t2)).^2);
	end
	d = find(squaredSum == min(squaredSum));
	d = d(1);
	delay = timeBins(d);

	if show,
		fig = figure;
		set(fig,'name','PETH Transition - Least Square Error','number','off');
		subplot(2,1,1);hold on;
		bar(timeBins,psth);
		plot([delay delay],ylim,'k','linestyle','--');
		title(['delay = ' num2str(delay) ' s']);
		ylabel('Occurrences');
		subplot(2,1,2);hold on;
		t1 = 1:d;t1 = t1';
		t2 = d+1:n;t2 = t2';
		plot(timeBins(t0),F0,'Color',[1 0 0],'Marker','none','LineStyle','-');
		plot(timeBins(t1+start-1),lambda1(d)*(t1)+F0(start),'Color',[0 0 0]);
		plot(timeBins(t2+start-1),lambda2(d)*(t2-d)+lambda1(d)*(d)+F0(start),'Color',[0 0 0]);
		plot([delay delay],ylim,'k','linestyle','--');
		xlabel('Time (s)');
		ylabel('Cumulative count');
	end

elseif strcmp(method,'pl')

	% Piecewise linear regression

	f = psth;

	n = length(f);
	t = 1:n;t = t';
	for theta = 2:n-2,
		t1 = 1:theta;t1 = t1';
		t2 = theta+1:n;t2 = t2';
		a1(theta) = LinearRegression(t1,f(t1));
		a2(theta) = LinearRegression(t2,f(t2),a1(theta)*theta);
		squaredSum(theta) = sum((f(t1)-a1(theta)*t1).^2) ...
			+ sum((f(t2)-a2(theta)*t2-a1(theta)*theta).^2);
	end
	squaredSum(1) = Inf;
	d = find(squaredSum == min(squaredSum));
	d = d(1);
	delay = timeBins(d);

	if show,
		fig = figure;
		set(fig,'name','PETH Transition - Piecewise Linear Fit','number','off');
		subplot(2,1,1);hold on;
		bar(timeBins,psth);
		plot([delay delay],ylim,'k','linestyle','--');
		title(['delay = ' num2str(delay) ' s']);
		ylabel('Occurrences');
		subplot(2,1,2);hold on;
		t1 = 1:d;t1 = t1';
		t2 = d+1:n;t2 = t2';
		plot(timeBins(t),f,'Color',[1 0 0],'Marker','o','LineStyle','none');
		plot(timeBins(t1),a1(d)*t1,'k');
		plot(timeBins(t2),a2(d)*t2+a1(d)*d,'k');
		plot([delay delay],ylim,'k','linestyle','--');
		xlabel('Time (s)');
		ylabel('Occurrences');
	end

end

function [a,b] = LinearRegression(x,y,b)

if nargin == 2, b = 0; end

y = y - b;
X = sum(x);
X2 = sum(x.^2);
Y = sum(y);
XY = sum(x.*y);
n = length(x);

b = (Y*X2-X*XY) / (n*X2-X.^2);
a = (n*XY-X*Y) / (n*X2-X.^2);
