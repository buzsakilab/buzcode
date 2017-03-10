function [consolidated,target] = ConsolidateIntervals(intervals,varargin)

%ConsolidateIntervals - Consolidate intervals.
%
% Consolidate overlapping intervals, e.g. replace [10,20] [15,25] with [10,25].
%
%  USAGE
%
%    [consolidated,target] = ConsolidateIntervals(intervals,<options>)
%
%    intervals      list of intervals
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'strict'      intervals with common bounds are consolidated ('off')
%                   or kept separate ('on') (default = 'off')
%     'epsilon'     intervals with close enough bounds (distance lesser than
%                   epsilon) are also consolidated (default = 0)
%    =========================================================================
%
%  OUTPUT
%
%    consolidated   consolidated intervals
%    target         for each original interval, the index of the consolidated
%                   interval to which it belongs (empty intervals yield NaN)
%
%  SEE
%
%    See also SubtractIntervals, ExcludeIntervals, InIntervals, Restrict,
%    FindInInterval, CountInIntervals, PlotIntervals.


% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
strict = 'off';
epsilon = 0;

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'strict',
			strict = lower(varargin{i+1});
			if ~isstring_FMAT(strict,'on','off'),
				error('Incorrect value for property ''strict'' (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
			end
		case 'epsilon',
			epsilon = varargin{i+1};
			if ~isdscalar(epsilon,'>0'),
				error('Incorrect value for property ''epsilon'' (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).']);
	end
end

original = intervals;

% Mark already consolidated intervals to avoid retesting them
done = logical(zeros(size(intervals(:,1))));

if strcmp(strict,'on'),
	for i = 1:length(intervals),
		if done(i), continue; end
		% Lower (L) and upper (U) interval bounds
		L = intervals(:,1);
		U = intervals(:,2);
		% Current interval is I = [l u], but we replace it with [l-e u+e] to take parameter 'epsilon' into account
		l = L(i)-epsilon;u = U(i)+epsilon;
		% Find all intervals that overlap with I:
		% 1) one of their bounds is strictly inside I
		% (their upper bound is greater than l, and their lower bound is lower than u)
		intersect = (U > l & L < u);
		% 2) they contain I
		if u == l,
 			% Special case: I is a singleton
 			intersect = intersect | (L < l & U > u);
		else
			intersect = intersect | (L <= l & U >= u);
		end
		% Determine smallest enclosing interval
		m = min(L(intersect));
		M = max(U(intersect));
		% Consolidate
		intervals(intersect,:) = repmat([m M],sum(intersect),1);
		done(intersect) = 1;
	end
else
	% (same as above, but replacing e.g. < with <=)
	for i = 1:length(intervals),
		if done(i), continue; end
		% Lower (L) and upper (U) interval bounds
		L = intervals(:,1);
		U = intervals(:,2);
		% Current interval is I = [l u], but we replace it with [l-e u+e] to take parameter 'epsilon' into account
		l = L(i)-epsilon;u = U(i)+epsilon;
		% Find all intervals that overlap with I:
		% 1) one of their bounds is inside I
		% (their upper bound is greater than or equal to l, and their lower bound is lower than or equal to u)
		intersect = (U >= l & L <= u);
		m = min(L(intersect));
		M = max(U(intersect));
		% Consolidate
		intervals(intersect,:) = repmat([m M],sum(intersect),1);
		done(intersect) = 1;
	end
end

% Sort intervals in ascending order (and store reordering information so we can reuse it later)
[intervals,order] = sortrows(intervals,1);

% Assign each consolidated interval an ID (in ascending order)
transitions = [1;find(diff(intervals(:,1))~=0)+1;length(intervals(:,1))];
for i = 1:length(transitions)-1,
	target(transitions(i):transitions(i+1)) = repmat(i,transitions(i+1)-transitions(i)+1,1);
end

%  Reorder consolidated interval IDs
target(order) = target;
target = target';

consolidated = unique(intervals,'rows');

% Remove empty intervals from output...
empty = diff(consolidated,1,2) < 0;
consolidated(empty,:) = [];
% ... and update target IDs
empty = diff(original,1,2) < 0;
[t,i] = sortrows([target empty]);
target(i) = t(:,1)-cumsum(t(:,2));

% Empty intervals belong to none
target(empty) = NaN;
