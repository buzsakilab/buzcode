%sca - Set current axes, without making the figure visible.
%
%  This function is useful in batches where figures must remain invisible. It
%  is equivalent to axes(h) but does not make the figure visible.
%
%  USAGE
%
%    sca(h)
%
%    h              axes handle
%
%  SEE ALSO
%
%    See also Hide, scf, StartBatch.
%

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function sca(h)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help sca">sca</a>'' for details).');
end

% Check parameters
if ~ishandle(h) || ~strcmp(get(h,'type'),'axes'),
  error('Incorrect axes handle (type ''help <a href="matlab:help sca">sca</a>'' for details).');
end

set(gcf,'CurrentAxes',h);
