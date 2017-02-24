function x = Clip(x,m,M)

%Clip - Clip values.
%
%  Clip values between specified extrema.
%
%  USAGE
%
%    clipped = Clip(x,min,max)
%
%    x              array to clip
%    min            minimum value
%    max            maximum value

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Clip">Clip</a>'' for details).');
end

x(x<m) = m;
x(x>M) = M;