%CountInIntervals - Count samples that fall in each of a list of intervals.
%
%  The equivalent Matlab code is trivial
%
%      for i = 1:length(intervals),
%        j = find(values>=intervals(i,1)&values<=intervals(i,2));
%        counts(i) = length(j);
%      end
%
%  but becomes extremely slow when dealing with very large lists.
%  This function can dramatically speed up things whenever one needs
%  to count values in a long list of intervals.
%
%  USAGE
%
%    counts = CountInIntervals(values,intervals)
%
%    values          values to test
%    intervals       [start,stop] pairs in ascending order
%
%  OUTPUT
%
%    counts         number of values that fall in each interval
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    InIntervals, Restrict, FindInInterval, PlotIntervals.

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
