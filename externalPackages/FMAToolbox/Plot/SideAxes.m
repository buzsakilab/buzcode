function h = SideAxes(a,location,s,varargin)

%SideAxes - Add side axes to existing axes.
%
%  USAGE
%
%    h = SideAxes(a,location,s,<options>)
%
%    Using cell arrays will overlay variable pairs.
%
%    a              optional target axes (default = gca)
%    s              size, as proportion of target axes
%    location       'top', 'bottom', 'left' or 'right'
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'gap'         gap between target and side axes, as proportion of total
%                   width or height (default = 0.1)
%    =========================================================================
%
%  OUTPUT
%
%    h              handle to the new axes


% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
gap = 0.1;

% Check number of parameters
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help SideAxes">SideAxes</a>'' for details).');
end

% Optional axes
if nargin == 2,
	s = location;
	location = a;
	a = gca;
elseif isstring_FMAT(a),
	varargin = {s,varargin{:}};
	s = location;
	location = a;
	a = gca;
end

% Check number of options
if mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help SideAxes">SideAxes</a>'' for details).');
end

% Check parameters
if ~ishandle(a),
	error('Incorrect axes (type ''help <a href="matlab:help SideAxes">SideAxes</a>'' for details).');
end
if ~isdscalar(s,'>0','<1'),
	error('Incorrect size (type ''help <a href="matlab:help SideAxes">SideAxes</a>'' for details).');
end
location = lower(location);
if ~isstring_FMAT(location,'top','bottom','left','right'),
	error('Incorrect location (type ''help <a href="matlab:help SideAxes">SideAxes</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SideAxes">SideAxes</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'gap',
			gap = varargin{i+1};
			if ~isdscalar(gap,'>=0','<1'),
				error('Incorrect value for property ''gap'' (type ''help <a href="matlab:help SideAxes">SideAxes</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SideAxes">SideAxes</a>'' for details).']);
	end
end

% Resize target axes to make room for side axes

main = a;
position = get(main,'position');
x = position(1);
y = position(2);
currentWidth = position(3); % current width
currentHeight = position(4); % current height

% The following variables will be used only when adding left or right side axes
horizontalGap = currentWidth * gap;
newWidth = currentWidth - horizontalGap;
newMainWidth = newWidth * (1-s);
newSideWidth = newWidth * s;

% The following variables will be used only when adding left or right side axes
verticalGap = currentHeight * gap;
newHeight = currentHeight - verticalGap;
newMainHeight = newHeight * (1-s);
newSideHeight = newHeight * s;

xLims = get(main,'xlim');
yLims = get(main,'ylim');

switch(location),
	case 'top',
		set(main,'position',[x y currentWidth newMainHeight]);
		h = axes('position',[x y+newMainHeight+verticalGap currentWidth newSideHeight]);
		xlim(h,xLims);
	case 'bottom',
		set(main,'position',[x y+newSideHeight+verticalGap currentWidth newMainHeight]);
		h = axes('position',[x y currentWidth newSideHeight]);
		xlim(h,xLims);
	case 'right',
		set(main,'position',[x y newMainWidth currentHeight]);
		h = axes('position',[x+newMainWidth+horizontalGap y newSideWidth currentHeight]);
		ylim(h,yLims);
	case 'left',
		set(main,'position',[x+newSideWidth+horizontalGap y newMainWidth currentHeight]);
		h = axes('position',[x y newSideWidth currentHeight]);
		ylim(h,yLims);
end
