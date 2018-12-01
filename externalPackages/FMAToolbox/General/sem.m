function s = sem(x,varargin)

%sem - Compute standard error of the mean (SEM).
%
%  USAGE
%
%    s = sem(x)
%
%    x              vector or matrix over which the SEM should be computed
%    dim            dimension to run SEM across (default = 1)

% Copyright (C) 2008-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help sem">sem</a>'' for details).');
end
if nargin == 1
    dim = 1;
else
    dim = varargin{1};
end
if ~isnumeric(x) & ~isdvector(x),
  error('Incorrect input - use vector or matrix (type ''help <a href="matlab:help sem">sem</a>'' for details).');
end

% if any(size(x)==1), x = x(:); end

% s = std(x,[],dim)./sqrt(size(x,1));
s = nanstd(x,[],dim) ./ sqrt(sum(~isnan(x),dim));