%clinspace - Linearly spaced vector of circular values (angles).
%
%  USAGE
%
%    x = clinspace(n,range)
%
%    n              number of values
%    range          optional:  1 for [-pi,pi[ (default)
%                              2 for [0,2pi[
%
%  SEE ALSO
%
%    See also isradians, wrap.
%

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function x = clinspace(n,range)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help clinspace">clinspace</a>'' for details).');
end

if nargin < 2,
	range = 1;
end

if ~isiscalar(n),
  error('Incorrect number of values (type ''help <a href="matlab:help clinspace">clinspace</a>'' for details).');
end

% Determine range
if range == 1,
	x = linspace(-pi,pi,n+1);
else
	x = linspace(0,2*pi,n+1);
end
x(end) = [];
