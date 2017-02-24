function line = UIAddLine(dir,varargin)

%UIAddLine - Interactively add horizontal/vertical line to existing plot.
%
%  USAGE
%
%    line = UIAddLine(dir,<options>)
%
%    dir            'v' for vertical, 'h' for horizontal
%    <options>      options for function <a href="matlab:help plot">plot</a>
%
%  NOTE
%
%    Click where the line should appear, then hit 'ENTER'.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help UIAddLine">UIAddLine</a>'' for details).');
end

hold on;
x = ginput;
if strcmp(lower(dir),'v'),
	line = plot([x(end,1) x(end,1)],ylim,varargin{:});
else
	line = plot(xlim,[x(end,2) x(end,2)],varargin{:});
end