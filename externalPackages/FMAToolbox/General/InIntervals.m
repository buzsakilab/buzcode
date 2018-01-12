function [status,interval,index] = InIntervals(values,intervals,varargin)

%InIntervals - Test which values fall in a list of intervals.
%
%  USAGE
%
%    [status,interval,index] = InIntervals(values,intervals,<options>)
%
%    values         values to test (these need not be ordered)
%    intervals      list of (start,stop) pairs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'verbose'     display information about ongoing processing
%                   (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    status         logical indices (1 = belongs to one of the intervals,
%                   0 = belongs to none)
%    interval       for each value, the index of the interval to which
%                   it belongs (0 = none)
%    index          for each value, its index in the interval to which
%                   it belongs (0 = none)
%
%  NOTE
%
%    If the intervals overlap, the outputs 'interval' and 'index' refer to the
%    last overlapping interval (i.e. if one value belongs to intervals #7 and #8,
%    it will be listed as belonging to interval #8).
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    Restrict, FindInInterval, CountInIntervals, PlotIntervals.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
verbose = false;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).');
end

% Check parameters
if ~isdmatrix(intervals) || size(intervals,2) ~= 2,
  error('Incorrect intervals (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).');
end

if isempty(values)
    warning('values is an empty vector, returning nothing..')
    status=[];
    interval=[];
    index=[];
    return
end

if size(values,1) == 1,
	values = values(:);
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'verbose',
			verbose = varargin{i+1};
			if ~isstring_FMAT(verbose,'on','off'),
				error('Incorrect value for property ''verbose'' (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).');
			end
			verbose = strcmp(verbose,'on');

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).']);
	end
end

[values,order] = sortrows(values(:,1));

% Max nb of digits for display
l = int2str(floor(log10(max(max(intervals*100))))+2);

% Determine if intervals overlap (in which case we must use a 'slow' algorithm)
i = sortrows(intervals,1);
di = i(2:end,1)-i(1:end-1,2);
overlap = any(di<0);
if ~overlap,
	% Fast algorithm: for the next interval, start from the end of the previous interval
	k = 2;
else
	% Slow algorithm: for the next interval, start from the beginning of the previous interval
	k = 1;
end

% Retrieve values in intervals
previous = 1;
n = size(intervals,1);
status = logical(zeros(size(values)));
interval = zeros(size(values));
index = zeros(size(values));
times = values;
for i = 1:n,
	from = intervals(i,1);
	to = intervals(i,2);
	timeString = sprintf(['%' l '.2f %' l '.2f (%' l '.2f)'],from,to,to-from);
	% Get values
	more = FindInInterval(values,[from to],previous);
	if ~isempty(more),
		previous = more(k); % See note above about algorithm
		nMore = more(2)-more(1)+1;
		interval(more(1):more(2)) = i;
		status(more(1):more(2)) = 1;
		index(more(1):more(2)) = (1:nMore);
	end
	if verbose, disp([timeString ' - ' int2str(nMore) ' values']); end
end

status(order) = status;
interval(order) = interval;
index(order) = index;
