%glinspace - Gamma-corrected linearly spaced vector.
%
%  USAGE
%
%    x = glinspace(d1,d2,gamma,n)
%
%    d1             start value
%    d2             stop value
%    gamma          gamma value
%    n              number of values
%
%  SEE ALSO
%
%    See also glinspace.
%

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function x = glinspace(d1,d2,gamma,n)

% Check number of parameters
if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help glinspace">glinspace</a>'' for details).');
end

% Check parameters
if ~isdscalar(d1),
  error('Incorrect start value (type ''help <a href="matlab:help glinspace">glinspace</a>'' for details).');
end
if ~isdscalar(d2),
  error('Incorrect stop value (type ''help <a href="matlab:help glinspace">glinspace</a>'' for details).');
end
if ~isdscalar(gamma) || gamma <= 0,
  error('Incorrect start value (type ''help <a href="matlab:help glinspace">glinspace</a>'' for details).');
end
if nargin < 3,
	n = 100;
end
if ~isiscalar(n),
  error('Incorrect number of values (type ''help <a href="matlab:help glinspace">glinspace</a>'' for details).');
end

if d1 > d2
	gamma = 1/gamma;
end

x = linspace(0,1,n) .^ (1/gamma);
x = x *(d2-d1) + d1;