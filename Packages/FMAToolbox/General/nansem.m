function s = nansem(x)

%nansem - Compute standard error of the mean (SEM), ignoring NaNs.
%
%  USAGE
%
%    s = nansem(x)
%
%    s              vector or matrix over which the sem should be computed

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help nansem">nansem</a>'' for details).');
end

if ~isdmatrix(x) & ~isdvector(x),
  error('Incorrect input - use vector or matrix (type ''help <a href="matlab:help nansem">nansem</a>'' for details).');
end

if any(size(x)==1), x = x(:); end

n = sum(~isnan(x));
s = nanstd(x)./sqrt(n);