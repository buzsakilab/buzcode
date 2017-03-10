%scf - Set current figure, without making it visible.
%
%  This function is useful in batches where figures must remain invisible. It
%  is equivalent to figure(h) but does not make h visible.
%
%  USAGE
%
%    scf(h)
%
%    h              figure handle
%
%  SEE ALSO
%
%    See also Hide, sca, StartBatch.
%

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function scf(h)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help scf">scf</a>'' for details).');
end

% Check parameters
if ~ishandle(h) || ~strcmp(get(h,'type'),'figure'),
  error('Incorrect figure handle (type ''help <a href="matlab:help scf">scf</a>'' for details).');
end

set(0,'CurrentFigure',h);
