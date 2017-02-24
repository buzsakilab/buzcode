function M = Array2Matrix(A)

%Array2Matrix - Transform an N-dimensional array into a matrix.
%
%  Each line in the output matrix lists the subscripts for all dimensions
%  (one per column) and the corresponding value in the array.
%

% Copyright (C) 2009-2011 by Michaël Zugaro
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
	error('Incorrect number of parameters (type ''help <a href="matlab:help Array">Array</a>2Matrix'' for details).');
end

d = length(size(A));

% Make a string of variable names to hold subscripts
%   e.g. for three dimensions: 'v1,v2,v3'
subscripts = '';
for i = 1:d,
	subscripts = [subscripts 'v' int2str(i) ','];
end
subscripts = [subscripts(1:end-1)];

% Make strings to transform indices to subscripts, and to create output matrix
%   e.g. for three dimensions: '[v1,v2,v3]=ind2sub(size(A),(1:length(A(:)))');'
%   and: 'M=[v1,v2,v3 A(:)];'
getSubscripts = ['[' subscripts ']=ind2sub(size(A),(1:length(A(:)))'');'];
setSubscripts = ['M=[' subscripts ' A(:)];'];

% Evaluate those strings
eval(getSubscripts);
eval(setSubscripts);
