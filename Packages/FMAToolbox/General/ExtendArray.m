function y = ExtendArray(x,s,fill)

%ExtendArray - Extend existing array to a given (larger) size.
%
%  USAGE
%
%    y = ExtendArray(x,s)
%
%    x              input array
%    s              output size, must have the same number of dimensions as x
%    fill           optional value for new elements (default: nan)
%

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check inputs
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ExtendArray">ExtendArray</a>'' for details).');
end
currentSize = size(x);
if length(currentSize) ~= length(s),
	error('Array and size do not have the same number of dimensions (type ''help <a href="matlab:help ExtendArray">ExtendArray</a>'' for details).');
end

% Create extended array, fill with requested value
if nargin < 3,
	y = nan(s);
else
	if ~isscalar(fill),
		error('Incorrect fill value (type ''help <a href="matlab:help ExtendArray">ExtendArray</a>'' for details).');
	end
	switch(fill),
		case 0,
			y = zeros(s);
		case 1,
			y = ones(s);
		case inf,
			y = inf(s);
		otherwise,
			y = fill*ones(s);
	end
end

% Copy old array into new array
indices = arrayfun(@(x) 1:x,currentSize,'uniformoutput',0);
y(indices{:}) = x;
