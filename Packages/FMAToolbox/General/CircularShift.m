function shifted = CircularShift(m,s)

%CircularShift - Shift matrix rows or columns circularly.
%
%  Shift each matrix row (or column) circularly by a different amount.
%
%  USAGE
%
%    shifted = CircularShift(m,s)
%
%    m              matrix to rotate
%    s              shift amount for each row (horizontal vector) or
%                   column (vertical vector)

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CircularShift">CircularShift</a>'' for details).');
end

[mm,nm] = size(m);
[ms,ns] = size(s);

% Check parameter sizes
if ~isivector(s),
	error('Second parameter is not a vector of integers (type ''help <a href="matlab:help CircularShift">CircularShift</a>'' for details).');
end
if mm ~= ms && nm ~= ns,
	error('Incompatible parameter sizes (type ''help <a href="matlab:help CircularShift">CircularShift</a>'' for details).');
end

% The algorithm below works along columns; transpose if necessary
s = -s(:)';
if ns == 1,
	m = m';
	[mm,nm] = size(m);
end

% Shift matrix S, where Sij is the vertical shift for element ij
shift = repmat(s,mm,1);

% Before we start, each element Mij has a linear index Aij. After circularly shifting the rows, it will have a linear index Bij.
% We now construct Bij.

% First, create matrix C where each item Cij = i (row number)
lines = repmat((1:mm)',1,nm);
% Next, update C so that Cij becomes the target row number (after circular shift)
lines = mod(lines+shift,mm);
lines(lines==0) = mm;
% Finally, transform Cij into a linear index, yielding Bij
indices = lines + repmat((0:nm-1).*mm,mm,1);

% Circular shift (reshape so that it is not transformed into a vector)
shifted = reshape(m(indices),mm,nm);

% Transpose back if necessary
if ns == 1,
	shifted = shifted';
end
