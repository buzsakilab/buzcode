function y = Insert(x,lines,where)

%Insert - Insert lines in a matrix.
%
%  USAGE
%
%    result = Insert(matrix,lines,indices)
%
%    array          matrix where the lines should be inserted
%    lines          list of values to insert
%    indices        list of (possibly repeated) matrix line numbers
%                   after which new lines should be inserted (use 0
%                   to insert before first line)
%
%  EXAMPLES
%
%    >> Insert([1;2;3;4;5],[-1;-2;-3],[0;3;5])
%
%    ans =
%
%      -1
%       1
%       2
%       3
%      -2
%       4
%       5
%      -3


% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

where = where(:);
mx = size(x,1);
nx = size(x,2);
ml = length(where);
if size(lines,1) == 1,
	lines = repmat(lines,ml,1);
end

% Lines in x must be moved to new indices in y
w = Accumulate(where+1,1,[mx+1,1]);
new = cumsum(w)+(1:mx+1)';
new(end) = [];

% Copy x into y
y = zeros(mx+ml,nx);
y(new,:) = x;

inserted = where+(1:ml)';
y(inserted,:) = lines;
