%FindInInterval - Find values that fall in a given interval.
%
%  The equivalent Matlab code is trivial
%
%      i = find(values(:,1)>=interval(1)&values(:,1)<=interval(2));
%      indices = [i(1),i(end)];
%
%  but becomes extremely slow when dealing with very large lists.
%  This function can dramatically speed up things whenever one needs
%  to repeatedly find values in a long list of intervals.
%
%  USAGE
%
%    indices = FindInInterval(values,interval,from)
%
%    values         values to test, sorted in ascending order
%    interval       [start,stop] pair
%    from           optional initial index (see example below)
%
%  OUTPUT
%
%    indices        indices of the first and last values that fall
%                   in the interval
%
%  EXAMPLE
%
%    This code assumes that the variable 'interval' contains a list
%    of non-overlapping intervals sorted in ascending order,
%    i.e. interval(i+1,1) >= interval(i,2).
%
%    previous = 1;
%    for i = 1:n,
%       % Find values within this window
%       j = FindInInterval(values,interval(i,:),previous);
%       previous = j(1);
%       % ... do whatever computations here ...
%    end
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    InIntervals, Restrict, CountInIntervals, PlotIntervals.

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
