function p = PlotTicks(xy,varargin)

%PlotTicks - Plot ticks.
%
%  USAGE
%
%    p = PlotTicks(x,y,<options>,<options2>)
%
%    x              list of tick abscissae
%    y              list of tick ordinates
%    <options>      optional list of property-value pairs (see table below)
%    <options2>     options for function <a href="matlab:help plot">plot</a>
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'direction'   direction: 'h' or 'v' (default = 'v')
%     'size'        tick size ([] or default = best guess)
%     'color'       either a single [R G B] color, or one per tick
%    =========================================================================
%
%  NOTE
%
%    When a single list of coordinates is provided, e.g. PlotTicks(t), these
%    are assumed to be abscissae or ordinates depending on 'direction', and
%    the ticks are plotted centered on 0.

% Copyright (C) 2008-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
 	error('Incorrect number of parameters (type ''help <a href="matlab:help PlotTicks">PlotTicks</a>'' for details).');
end

% Defaults
direction = 'v';
yx = 0;
color = [];
s = [];
v = {};

% Does parameter list include both x and y?
both = false;
if ~isempty(varargin) && isnumeric(varargin{1}),
	x = xy;
	y = varargin{1};
	varargin = {varargin{2:end}};
	both = true;
end

% Parse parameter list
n = 1;
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotTicks">PlotTicks</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'size',
			s = varargin{i+1};
			if ~isdscalar(s) && ~isdvector(s),
				error('Incorrect value for property ''size'' (type ''help <a href="matlab:help PlotTicks">PlotTicks</a>'' for details).');
			end

		case 'direction',
			direction = lower(varargin{i+1});
			if ~isstring_FMAT(direction,'h','v'),
				error('Incorrect value for property ''direction'' (type ''help <a href="matlab:help PlotTicks">PlotTicks</a>'' for details).');
			end

		case 'color',
			color = varargin{i+1};
			if ~isdmatrix(color,'>=0','<=1') || size(color,2) ~= 3,
				error('Incorrect value for property ''color'' (type ''help <a href="matlab:help PlotTicks">PlotTicks</a>'' for details).');
			end

		otherwise,
			v = {varargin{i:end}};
			if ~isa(v,'cell'), v = {v}; end
			break;
	end
end

if ~both,
	if strcmp(direction,'v'),
		x = xy;
		y = yx;
	else
		x = yx;
		y = xy;
	end
end

% Check list lengths
if ~isdvector(x) | ~isdvector(y) | (~isdscalar(x) & ~isdscalar(y) & any(size(x) ~= size(y))),
 	error('Incorrect abscissae/ordinates (type ''help <a href="matlab:help PlotTicks">PlotTicks</a>'' for details).');
else
	if isscalar(x),
		x = x*ones(size(y));
	else
		x = x(:);
	end
	if isscalar(y),
		y = y*ones(size(x));
	else
		y = y(:);
	end
end

% Default size
if isempty(s),
	if strcmp(direction,'v'),
		s = min(diff(unique(y)))*0.75;
	else
		s = min(diff(unique(x)))*0.75;
	end
	if isempty(s), s = 0.75; end
end

s = s(:);
if length(s) == 1,
	s = repmat(s,length(x),1);
elseif any(size(x) ~= size(s)),
 	error('Inconsistent numbers of abscissae/ordinates and sizes (type ''help <a href="matlab:help PlotTicks">PlotTicks</a>'' for details).');
end

hold on;
if strcmp(direction,'v'),
	if isempty(color),
		for i = 1:length(x),
			plot([x(i) x(i)],[y(i)-s(i)/2 y(i)+s(i)/2],v{:});
		end
	else
		if size(color,1) == 1, color = repmat(color,length(x),1); end
		for i = 1:length(x),
			plot([x(i) x(i)],[y(i)-s(i)/2 y(i)+s(i)/2],'color',color(i,:),v{:});
		end
	end
else
	if isempty(color),
		for i = 1:length(x),
			plot([x(i)-s(i)/2 x(i)+s(i)/2],[y(i) y(i)],v{:});
		end
	else
		if size(color,1) == 1, color = repmat(color,length(x),1); end
		for i = 1:length(x),
			plot([x(i)-s(i)/2 x(i)+s(i)/2],[y(i) y(i)],'color',color(i,:),v{:});
		end
	end
end
