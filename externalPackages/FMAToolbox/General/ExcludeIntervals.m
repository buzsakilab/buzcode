function [intervals,indices] = ExcludeIntervals(intervals1,intervals2,varargin)

%ExcludeIntervals - Exclude intersecting intervals.
%
% Given two lists of  intervals, exclude from the first list those intervals
% that intersect with one or more of the intervals in the second list.
%
%  USAGE
%
%    [intervals,indices] = ExcludeIntervals(intervals1,intervals2,<options>)
%
%    intervals1     list of reference intervals
%    intervals2     list of intervals to exclude
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'strict'      intervals with common bounds are as intersecting ('on')
%                   or disjoint ('off') (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    intervals      non-intersecting intervals
%    indices        indices of these intervals in the original list
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, InIntervals, Restrict,
%    FindInInterval, CountInIntervals, PlotIntervals.


% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
strict = 'off';

if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ExcludeIntervals">ExcludeIntervals</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ExcludeIntervals">ExcludeIntervals</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help <a href="matlab:help ExcludeIntervals">ExcludeIntervals</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'strict',
			strict = lower(varargin{i+1});
			if ~isstring_FMAT(strict,'on','off'),
				error('Incorrect value for property ''strict'' (type ''help <a href="matlab:help ExcludeIntervals">ExcludeIntervals</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ExcludeIntervals">ExcludeIntervals</a>'' for details).']);
	end
end

intervals = [];
indices = [];
if strcmp(strict,'on'),
	for i = 1:length(intervals1),
		% List all intervals (to exclude) that overlap with the current interval
		intersect = (intervals1(i,1) <= intervals2(:,1) & intervals1(i,2) > intervals2(:,1)) ...
			| (intervals1(i,1) >= intervals2(:,1) & intervals1(i,2) <= intervals2(:,2)) ...
			| (intervals1(i,1) < intervals2(:,2) & intervals1(i,2) >= intervals2(:,2));
		% Include current interval?
		if sum(intersect) == 0,
			intervals = [intervals;intervals1(i,:)];
			indices = [indices;i];
		end
	end
else
	for i = 1:length(intervals1),
		% List all intervals (to exclude) that overlap with the current interval
		intersect = (intervals1(i,1) <= intervals2(:,1) & intervals1(i,2) >= intervals2(:,1)) ...
			| (intervals1(i,1) >= intervals2(:,1) & intervals1(i,2) <= intervals2(:,2)) ...
			| (intervals1(i,1) <= intervals2(:,2) & intervals1(i,2) >= intervals2(:,2));
		% Include current interval?
		if sum(intersect) == 0,
			intervals = [intervals;intervals1(i,:)];
			indices = [indices;i];
		end
	end
end
