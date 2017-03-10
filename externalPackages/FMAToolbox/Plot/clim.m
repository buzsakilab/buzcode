function c = clim(handle,clims)

%clim - Get or set color scaling limits for current axes (or all figure axes).
%
%  USAGE
%
%    clims = clim
%    clim(h)
%    clim([c0 c1])
%    clim(fig,[c0 c1])
%
%    h              optional axes or figure handle (default = gca)
%    c0, c1         optional scaling minimum and maximum


% Copyright (C) 2009-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

switch(nargin),
	case 0,
		% Return clim for current axes
		c = get(gca,'CLim');
	case 1,
		if isdscalar(handle),
			% Return clim for the given axes
			if ~ishandle(handle),
				error('Incorrect handle (type ''help <a href="matlab:help clim">clim</a>'' for details).');
			end
			c = get(handle,'CLim');
		else
			% Set clim for current axes
			clims = handle;
			if ~isdvector(clims,'<'),
				error('Incorrect color limits (type ''help <a href="matlab:help clim">clim</a>'' for details).');
			end
			set(gca,'CLim',clims);
		end
	case 2,
		if ~isdscalar(handle) || ~ishandle(handle),
			error('Incorrect handle (type ''help <a href="matlab:help clim">clim</a>'' for details).');
		end
		if ~isdvector(clims,'<'),
			error('Incorrect color limits (type ''help <a href="matlab:help clim">clim</a>'' for details).');
		end
		% Set clim for axes or figure
		c = clims;
		if strcmp(get(handle,'type'),'figure'),
			% Set clim for all axes of a figure
			children = get(handle,'children');
			for i = 1:length(children),
				child = children(i);
				if strcmp(get(child,'type'),'axes'),
					set(child,'CLim',clims);
				end
			end
		else
			% Set clim for these axes
			set(handle,'CLim',clims);
		end
	end
end