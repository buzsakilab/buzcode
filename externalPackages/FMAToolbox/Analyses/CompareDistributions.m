function [h,stats] = CompareDistributions(group1,group2,varargin)

%CompareDistributions - Compare two N-dimensional distributions (e.g. prospective place fields).
%
%  This non-parametric test determines in which dimensions two distributions are
%  significantly different. For example, the two distributions could be
%  leftward vs rightward trajectories in a T maze (each 'dimension' of these
%  distributions corresponds to a spatial bin), and the test would determine
%  where the trajectories are significantly different (diverge).
%
%  This test is based on the bootstrap method developed in Fujisawa et al.
%  (2008). This compares the observed difference between the two distribution
%  averages vs the averages and confidence intervals for surrogate data (where
%  individual observations are shuffled across groups). Both pointwise and
%  global confidence intervals are computed.
%
%  USAGE
%
%    [h,stats] = CompareDistributions(group1,group2,<options>)
%
%    group1,group2  values for the two groups: lines are observations, columns
%                   are dimensions (e.g. spatial bins)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nShuffles'   number of shuffles to estimate the distribution
%                   (default = 5000)
%     'alpha'       confidence level (default = 0.05)
%     'max'         maximum number of iterations for global confidence
%                   intervals (default = 6)
%     'tolerance'   maximum difference (in %) between target and computed
%                   global alpha (default = 0.8)
%     'tail'        one or two-tailed ditributions (default = two)
%     'show'        show figure (default 'off')
%     'verbose'     show information about ongoing processing (default 'off')
%    =========================================================================
%
%  OUTPUT
%
%    h                  h = 1 if the null hypothesis can be rejected
%    stats.observed     observed difference between the means of the two groups
%    stats.null         expected difference between the means of the two groups
%                       under the null hypothesis
%    stats.pointwise    pointwise confidence intervals at alpha level
%    stats.global       global confidence intervals at alpha level
%    stats.above        logical vector indicating where the null hypothesis can
%                       be rejected because the observed difference exceeds the
%                       upper confidence limits
%    stats.below        logical vector indicating where the null hypothesis can
%                       be rejected because the observed difference lies below
%                       the lower confidence limits
%    stats.alpha        successive pointwise alpha values used to target
%                       the required global alpha level (one for each iteration)
%    stats.p            successive resulting global alpha values (one for each iteration)
%

% Copyright (C) 2010-2011 by Erika Cerasti, 2013-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
nShuffles = 5000;
alpha = 0.05;
tail = 'two';
show = false;
maxIterations = 6;
tolerance = 0.8;
verbose = false;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
    error('Incorrect number of parameters (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
end

if isempty(group1) || isempty(group2), return; end

% Check parameter sizes
if size(group1,2) ~= size(group2,2),
    error('The two groups do not have the same number of columns (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
end

isNDimensional = size(group1,2) > 1; % N dimensions (N>1) => compute global confidence intervals

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'nshuffles',
			nShuffles = lower(varargin{i+1});
			if ~isdscalar(nShuffles,'>0'),
				error('Incorrect value for property ''nShuffles'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'alpha',
			alpha = varargin{i+1};
			if ~isdscalar(alpha,'>0','<1'),
				error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'max',
			maxIterations = varargin{i+1};
			if ~isdscalar(maxIterations,'>1'),
				error('Incorrect value for property ''max'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'tolerance',
			tolerance = varargin{i+1};
			if ~isdscalar(tolerance,'>0'),
				error('Incorrect value for property ''tolerance'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'tail',
			tail = varargin{i+1};
			if ~isstring(tail,'one','two'),
				error('Incorrect value for property ''tail'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
			show = strcmp(show,'on');
		case 'verbose',
			verbose = varargin{i+1};
			if ~isstring(verbose,'on','off'),
				error('Incorrect value for property ''verbose'' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).');
			end
			verbose = strcmp(verbose,'on');
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CompareDistributions">CompareDistributions</a>'' for details).']);
	end
end

if verbose,
	disp('Shuffling the data...');
end

n1 = size(group1,1);
n2 = size(group2,1);

differences = nan(nShuffles,size(group1,2));

g = [group1;group2];
for i = 1:nShuffles,
	
   % Shuffle the data
	shuffled = randperm(n1+n2);
	shuffled1 = g(shuffled(1:n1),:);
	shuffled2 = g(shuffled(n1+1:end),:);

	% Compute the mean difference
	m1 = mean(shuffled1);
	m2 = mean(shuffled2);
	differences(i,:) = m1 - m2;

	if verbose,
		step = i - floor(i/1000)*1000;
		if step == 0,
			disp([' Iteration # ' int2str(i) '/' int2str(nShuffles)]);
		end
	end

end

if verbose,
	disp('Computing confidence intervals...');
end

nBins = size(group1,2);
deviation = 1;
iteration = 0;
while deviation*100 > tolerance,

	% Check number of iterations
	iteration = iteration + 1;
	if iteration > maxIterations,
		warning(['Reached maximum number of iterations (' int2str(maxIterations) ')']);
		break;
	end

	% Update quantiles based on the value of alpha computed during the previous iteration
	if strcmp(tail,'one'),
		quantiles = 1-alpha;
	else
		quantiles = [alpha/2 1-(alpha/2)];
	end

	% Pointwise confidence intervals at alpha level, using the value of alpha computed during the previous iteration
	confidenceIntervals = quantile(differences,quantiles);
	confidenceIntervals = confidenceIntervals([2 1],:);

	if ~isNDimensional,
		pointwise = confidenceIntervals;
		break;
	end
	
	% Global confidence intervals, using the pointwise intervals
	significant = differences > repmat(confidenceIntervals(1,:),nShuffles,1) | differences < repmat(confidenceIntervals(2,:),nShuffles,1);
	n = sum(any(significant,2));
	p = n/nShuffles;
	pGlobal(iteration) = p;
	alphaGlobal(iteration) = alpha;

	% Update alpha for next iteration
	deviation = abs(p-alpha);
	if iteration == 1,
		pointwise = confidenceIntervals;
		alpha = 0.003;
	else
		s = sign(p-alpha);
		if deviation > 0.05,
			alpha = alpha - s*deviation*0.00075;
		else
			alpha = alpha - s*deviation*0.033;
		end
		if alpha < 0.001,
			warning('alpha is becoming too small, global confidence bands will not be accurate');
			break;
		end
	end

end

% Statistics
stats.observed = mean(group1) - mean(group2);
stats.null = mean(differences,1);
stats.pointwise = pointwise;
stats.global = [];
stats.alpha = [];
stats.p = [];
if(iteration > 1)
    stats.global = confidenceIntervals;
    stats.alpha = alphaGlobal;
    stats.p = pGlobal;
end

% Determine significant segments

if ~isNDimensional,

    h = stats.observed > stats.pointwise(1) | stats.observed < stats.pointwise(2);
    stats.above = [];
    stats.below = [];

else

	% Where does the observed difference exceed the pointwise upper confidence limit? (etc.)
	abovePointwise = stats.observed > stats.pointwise(1,:);
	belowPointwise = stats.observed < stats.pointwise(2,:);
	aboveGlobal = stats.observed > stats.global(1,:);
	belowGlobal = stats.observed < stats.global(2,:);
	above = abovePointwise & aboveGlobal;
	below = belowPointwise & belowGlobal;

	% Segments above the upper limit: a segment is considered significant if exceeds the global upper limit,
	% but its extent is determined using the pointwise upper limit

	if any(above),
		% Points where global upper limit is exceeded (when several points are adjacent, keep only the first one)
		start = find(diff([0 above])==1);
		for i = start,
			% Examine points to the left: do they exceed the pointwise upper limit?
			j = 1;
			while i-j ~= 0 && abovePointwise(i-j),
				above(i-j) = 1;
				j = j+1;
			end
			% Examine points to the right: do they exceed the pointwise upper limit?
			j = 1;
			while i+j <= length(abovePointwise) && abovePointwise(i+j),
				above(i+j) = 1;
				j = j+1;
			end
		end
	end
	stats.above = above;

	% Segments below the lower limit
	if any(below),
		start = find(diff([0 below])==1);
		for t=1: length(start)
			i=start(t);
			j=1;
			while(i-j)~=0 && (belowPointwise(i-j)==1)
					below(i-j)=1;
					j=j+1;
			end
			j=1;
			while(i+j)<=length(belowPointwise) && (belowPointwise(i+j)==1)
					below(i+j)=1;
					j=j+1;
			end
		end
	end
	stats.below = below;

	% Reject H0?
	h = any(above) | any(below);

end

if isNDimensional && show,

	orange = [1 0.4 0];
	lightGreen = [0.6 0.9 0.2];
	darkGreen = [0.15 0.35 0.15];
	
	x = 1:length(stats.observed);
	figure; hold on;
	plot(x,stats.global,'o-','MarkerSize',4,'Color',darkGreen,'MarkerFaceColor',darkGreen,'MarkerEdgeColor',darkGreen);
	plot(x,stats.pointwise,'o-','MarkerSize',4,'Color',lightGreen,'MarkerFaceColor',lightGreen,'MarkerEdgeColor',lightGreen);
	plot(x,stats.observed,'Color',orange,'LineWidth',2);
	PlotHVLines(0,'h','--','Color','k');
	ylabel('Mean Differences');
	xlabel('Bins');
	
	PlotIntervals(ToIntervals(stats.above|stats.below),'rectangles');

end
