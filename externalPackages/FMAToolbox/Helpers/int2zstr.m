%int2zstr - Convert integer to zero-padded string.
%
%  USAGE
%
%    s = int2zstr(n,digits)
%
%    n              integer to convert
%    digits         optional total number of digits
%

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function s = int2zstr(n,digits)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help int2zstr">wrap</a>'' for details).');
end

% Convert to string
s = int2str(n);

% Pad if necessary
if nargin < 2, return; end
l = digits-length(s);
if l <= 0, return; end
s = [ones(1,l)*'0' s];
