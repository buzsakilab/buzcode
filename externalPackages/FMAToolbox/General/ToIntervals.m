function intervals = ToIntervals(x,in)

%ToIntervals - Convert logical vector to a list of intervals.
%
%  USAGE
%
%    intervals = ToIntervals(x,in)
%
%    x              values, e.g. timestamps
%    in             logical vector (1 = x value is inside, 0 = x value is outside)
%
%  NOTE
%
%    Values can also be omitted, in which case the intervals are defined in terms
%    of indices in the logical vector (see Example below).
%
%  EXAMPLES
%
%
%    x = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1];
%    in = [0 0 1 1 1 0 0 0 1 1 0];
%    ToIntervals(x,in)
%
%     ans =
%               0.3      0.5
%               0.9        1
%
%    ToIntervals(in)
%
%     ans =
%               3      5
%               9     10
%
%  SEE
%
%    See also ConsolidateIntervals, ExcludeIntervals, InIntervals, Restrict,
%    FindInInterval, CountInIntervals, PlotIntervals.
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ToIntervals">ToIntervals</a>'' for details).');
end

if nargin == 1,
	in = x;
end

if ~islvector(in),
	error('Incorrect logical vector (type ''help <a href="matlab:help ToIntervals">ToIntervals</a>'' for details).');
end

if in(end) == 1, in(end+1) = 0; end
in = in(:);
din = diff([0;in]);

start = din == 1;
stop = din == -1;

intervals = [find(start) find(stop)-1];

if nargin >=2,
	if ~isdvector(x),
		error('Incorrect x values (type ''help <a href="matlab:help ToIntervals">ToIntervals</a>'' for details).');
	end
	intervals = x(intervals);
end
