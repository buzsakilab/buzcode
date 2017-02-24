function M = Array2PagedMatrix(A)

%Array2PagedMatrix - Transform an N-dimensional array into a paged matrix.
%
%  Vertically concatenate all pages in the array, and add N-1
%  columns listing the subscripts for all dimensions except the
%  second dimension (which is implicit, since it becomes the
%  second dimension of the output matrix).
%
%  Consider the following example. Let A be the array:
%
%    A(:,:,1,1) =
%
%       1     2     3     4     5
%       6     7     8     9     0
%
%
%    A(:,:,2,1) =
%
%       0     9     8     7     6
%       5     4     3     2     1
%
%
%    A(:,:,1,2) =
%
%       4     4     4     4     4
%       5     5     5     5     5
%
%
%    A(:,:,2,2) =
%
%       6     6     6     6     6
%       7     7     7     7     7
%
%
%    A(:,:,1,3) =
%
%       8     8     8     8     8
%       9     9     9     9     9
%
%
%    A(:,:,2,3) =
%
%       0     0     0     0     0
%       6     5     6     5     6
%
%  Array2PagedMatrix concatenates each of the above pages vertically, and adds three
%  columns listing the first, third and fourth dimensions, yielding:
%
%    B =
%
%       1     2     3     4     5     1     1     1
%       6     7     8     9     0     2     1     1
%       0     9     8     7     6     1     2     1
%       5     4     3     2     1     2     2     1
%       4     4     4     4     4     1     1     2
%       5     5     5     5     5     2     1     2
%       6     6     6     6     6     1     2     2
%       7     7     7     7     7     2     2     2
%       8     8     8     8     8     1     1     3
%       9     9     9     9     9     2     1     3
%       0     0     0     0     0     1     2     3
%       6     5     6     5     6     2     2     3


% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%  A(:,:,1,1) = [1 2 3 4 5;6 7 8 9 0];
%  A(:,:,2,1) = [0 9 8 7 6;5 4 3 2 1];
%  A(:,:,1,2) = [4 4 4 4 4;5 5 5 5 5];
%  A(:,:,2,2) = [6 6 6 6 6;7 7 7 7 7];
%  A(:,:,1,3) = [8 8 8 8 8;9 9 9 9 9];
%  A(:,:,2,3) = [0 0 0 0 0;6 5 6 5 6];

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Array">Array</a>2PagedMatrix'' for details).');
end

if length(size(A)) <= 2
	M = A;
	return
end

% The array A contains k pages, each of which is a m x n matrix. Let d be the number of
% dimensions. The output matrix should have size k*m x (n+d-1)
m = size(A,1);
n = size(A,2);
d = length(size(A));
k = prod(size(A))/(m*n);
M = zeros((k*m),(n+d-1));

% This string will be evaluated later to determine for each page the subscripts for the d extra dimensions
% (we need to use this trick because we do not know how many extra dimensions there are)
getSubscripts = '[';
for i = 1:d,
	getSubscripts = [getSubscripts 'v' int2str(i) ','];
end
getSubscripts = [getSubscripts(1:end-1) ']=ind2sub(size(A),i*s);'];

% This string will be evaluated later to set the last d columns (corresponding to the subscripts for the d extra dimensions)
% (we need to use this trick because we do not know how many extra dimensions there are)
setSubscripts = 'M((i-1)*m+1:i*m,n+2:end)=repmat([';
for i = 3:d,
	setSubscripts = [setSubscripts 'v' int2str(i) ','];
end
setSubscripts = [setSubscripts(1:end-1) '],m,1);'];

s = m*n;
for i = 1:k,
	% Get i-th page
	M((i-1)*m+1:i*m,1:n) = reshape(A((i-1)*s+1:i*s),m,n);
	% Set the last d columns (corresponding to the subscripts for the extra dimensions)
	eval(getSubscripts);
	eval(setSubscripts);
end

M(:,n+1) = repmat((1:m)',k,1);