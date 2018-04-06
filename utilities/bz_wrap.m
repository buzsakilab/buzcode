%wrap - Set radian angles in range [0,2pi] or [-pi,pi].
%
%  USAGE
%
%    y = wrap(x,offset)
%
%    x              angles in radians
%    offset         phase angle [0-2*pi] that circularly shifts output data
%    
%
%  SEE ALSO
%
%    See also isradians, clinspace.
%

% Copyright (C) 2010-2011 by Michaël Zugaro
% edited by David Tingley, 2017 to add phase offset as input
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function y = bz_wrap(x,offset)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help wrap">wrap</a>'' for details).');
end

% if nargin < 2,
% 	range = 1;
% end

if ~isa(x,'double'), y = []; return; end

% Determine angle in [0,2*pi]
y = mod(mod(x,2*pi)+offset,2*pi);

% Change range if necessary
% if range == 1,
% 	change = y > pi;
% 	y(change) = y(change)-2*pi + offset;
% end