function c = clim(clims)

%clim - Get or set color scaling limits for current axes.
%
%  USAGE
%
%    clims = clim
%    clim([c0 c1])
%
%    c0, c1         scaling minimum and maximum


% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin ==1,
	set(gca,'CLim',clims);
end
c = get(gca,'CLim');

