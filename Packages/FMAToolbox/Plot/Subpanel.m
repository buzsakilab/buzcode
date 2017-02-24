function panel = Subpanel(m,n,p,varargin)

%Subpanel - Add an (invisible) panel container to a figure.
%
%  USAGE
%
%    h = Subpanel(m,n,p,<options>)
%
%    m,n,p          panel position as in <a href="matlab:help subplot">subplot</a>
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'parent'      panel parent handle (default = gcf)
%    =========================================================================
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
parent = [];

% Check number of parameters
if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Subpanel">Subpanel</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Subpanel">Subpanel</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'parent',
			parent = varargin{i+1};
			if ~ishandle(parent),
				error('Incorrect value for property ''parent'' (type ''help <a href="matlab:help Subpanel">Subpanel</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Subpanel">Subpanel</a>'' for details).']);
	end
end

if isempty(parent), parent = gcf; end

% Create panel (where the corresponding subplot would be)
position = [0 0 1 1];
x0 = position(1);
dx = position(3)/n;
y0 = position(2);
dy = position(4)/m;
p = p - 1;
x = rem(p,n);
y = floor(p/n);
X1 = min(x);
Y1 = min(y);
X2 = max(x);
Y2 = max(y);
P = [x0+X1*dx y0+(m-Y2-1)*dy (X2-X1+1)*dx (Y2-Y1+1)*dy];
panel = uipanel('position',P,'parent',parent);

% Change appearance
if strcmp(get(parent,'type'),'figure'),
	color = get(parent,'color');
elseif strcmp(get(parent,'type'),'uipanel'),
	color = get(parent,'backgroundcolor');
else
	color = [1 1 1]*0.8;
end
set(panel,'bordertype','none','backgroundcolor',color);
