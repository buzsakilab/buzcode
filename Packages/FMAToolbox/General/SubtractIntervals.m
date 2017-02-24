function [intervals,indices] = SubtractIntervals(intervals1,intervals2,varargin)

%SubtractIntervals - Subtract intervals.
%
% Given two lists of intervals, subtract from each interval in the first
% list its intersection with each of the intervals in the second list.
%
%  USAGE
%
%    [intervals,indices] = SubtractIntervals(intervals1,intervals2,<options>)
%
%    intervals1     list of reference intervals
%    intervals2     list of intervals to subtract
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
%    intervals      resulting intervals
%    indices        indices of these intervals in the original list
%
%  SEE
%
%    See also ConsolidateIntervals, ExcludeIntervals, InIntervals, Restrict,
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
	error('Incorrect number of parameters (type ''help <a href="matlab:help SubtractIntervals">SubtractIntervals</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help SubtractIntervals">SubtractIntervals</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help <a href="matlab:help SubtractIntervals">SubtractIntervals</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'strict',
			strict = lower(varargin{i+1});
			if ~isstring_FMAT(strict,'on','off'),
				error('Incorrect value for property ''strict'' (type ''help <a href="matlab:help SubtractIntervals">SubtractIntervals</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SubtractIntervals">SubtractIntervals</a>'' for details).']);
	end
end

indices = (1:size(intervals1,1))';
for i = 1:size(intervals2,1),
	% Trim all intervals in intervals1 that left-intersect with intervals2(i)
	intersect = intervals2(i,1) <= intervals1(:,1) & intervals2(i,2) >= intervals1(:,1);
	if sum(intersect) ~= 0, intervals1(intersect,1) = intervals2(i,2); end
	% Trim all intervals in intervals1 that right-intersect with intervals2(i)
	intersect = intervals2(i,1) <= intervals1(:,2) & intervals2(i,2) >= intervals1(:,2);
	if sum(intersect) ~= 0, intervals1(intersect,2) = intervals2(i,1); end
	% Remove all intervals in intervals1 that are contained in intervals2(i)
	inside = intervals2(i,1) <= intervals1(:,1) & intervals2(i,2) >= intervals1(:,2);
	if sum(inside) ~= 0, intervals1(inside,:) = []; indices(inside,:) = []; end
	% Split all intervals in intervals1 that contain intervals2(i)
	contain = intervals2(i,1) >= intervals1(:,1) & intervals2(i,2) <= intervals1(:,2);
	n = sum(contain);
	if n ~= 0,
		first = [intervals1(contain,1) repmat(intervals2(i,1),n,1)];
		second = [repmat(intervals2(i,2),n,1) intervals1(contain,2)];
		intervals1(contain,:) = first;
		intervals1 = Insert(intervals1,second,find(contain));
		indices = Insert(indices,indices(contain),find(contain));
	end
end
intervals = intervals1;
