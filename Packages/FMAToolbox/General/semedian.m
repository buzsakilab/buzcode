function s = semedian(x,varargin)

%semedian - Compute standard error of the median.
%
%  USAGE
%
%    s = semedian(x)
%
%    x              vector or matrix over which the error should be computed

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help semedian">semedian</a>'' for details).');
end
if ~isdmatrix(x) & ~isdvector(x),
  error('Incorrect input - use vector or matrix (type ''help <a href="matlab:help semedian">semedian</a>'' for details).');
end

if any(size(x)==1), x = x(:); end

n = size(x,1);
m = repmat(median(x),n,1);
s = sqrt( sum((x-m).^2) / (n*(n-1)) );
