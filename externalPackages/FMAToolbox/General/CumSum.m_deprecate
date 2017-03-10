function s = CumSum(data,stops)

%CumSum - Cumulative sum of elements. Partial sums can also be computed.
%
%  USAGE
%
%    sum = CumSum(data,stops)
%
%    data           data to sum
%    stops          optional logical indices where sum should be restarted
%
%  EXAMPLE
%
%    % Simple cumulative sum
%    s = CumSum([1 4 6 2 13]);
%
%    % Partial cumulative sums
%    s = CumSum([1 4 6 2 13 2 4 6 5 10 1],[1 0 0 0 0 1 0 0 0 0 0])
%

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help CumSum">CumSum</a>'' for details).');
end
if ~isvector(data),
	error('Parameter ''data'' is not a vector (type ''help <a href="matlab:help CumSum">CumSum</a>'' for details).');
end
data = data(:);
if nargin == 2,
	if ~isvector(stops),
		error('Parameter ''stops'' is not a vector (type ''help <a href="matlab:help CumSum">CumSum</a>'' for details).');
	end
	stops = logical(stops(:));
	if length(stops) ~= length(data),
		error('Parameters ''data'' and ''stops'' have different lengths (type ''help <a href="matlab:help CumSum">CumSum</a>'' for details).');
	end
end

s = cumsum(data)
if nargin == 1, return; end

data(stops) = 2*data(stops)-s(stops);
s = cumsum(data);
