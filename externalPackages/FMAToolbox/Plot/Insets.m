function [a,varargout] = Insets(location,s,varargin)

%Insets - Create insets in current axes.
%
%  Create insets (= small figures) in current axes. Insets are laid out side
%  by side from left to right. The number of insets is determined by the
%  number of output parameters.
%
%  USAGE
%
%    [a,b,...] = Insets(location,size,<options>)
%
%    location       'topleft', 'top', 'topright', 'left', 'center', 'right',
%                   'bottomleft', 'bottom', or 'bottomright'
%    size           total [h v] size of insets, as proportions of current axes
%                   (after subtracting margins)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'margin'      margin around insets, as a proportion of figure size
%                   (default = 0.1)
%    =========================================================================
%
%  EXAMPLE
%
%    % Create a figure and add the main plot
%    figure;plot(x,y);
%    % Create 3 insets in the top right corner, leaving 10% of the available
%    % width and height for margins (default). In total, the insets will occupy
%    % half of the remaining width, and one quarter of the remaining height.
%    [a,b,c] = Insets('topright',[0.5 0.25]);
%    % Plot data in each of the insets
%    axes(a);plot(x,z1);
%    axes(b);plot(x,z2);
%    axes(c);plot(x,z3);

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
margin = 0.1;

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Insets">Insets</a>'' for details).');
end

% Check parameters
if any(s>1|s<0),
  error('Sizes must be expressed as proportion [0..1] of total figure size (type ''help <a href="matlab:help Insets">Insets</a>'' for details).');
end
location = lower(location);
if ~isstring_FMAT(location,'topleft','top','topright','left','center','right','bottomleft','bottom','bottomright'),
  error('Incorrect location (type ''help <a href="matlab:help Insets">Insets</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Insets">Insets</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'margin',
			margin = lower(varargin{i+1});
			if ~isdscalar(margin,'>=0','<=1'),
				error('Incorrect value for property ''margin'' (type ''help <a href="matlab:help Insets">Insets</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Insets">Insets</a>'' for details).']);
	end
end


% Position and size of current axes
p = get(gca,'position');
axesBottom = p(2);
axesTop = p(2)+p(4);
axesLeft = p(1);
axesRight = p(1)+p(3);
axesWidth = p(3);
axesHeight = p(4);

% Inset margins
hMargin = margin*axesWidth;
vMargin = margin*axesHeight;

% Size of each mini figure
n = nargout;
dx = s(1)/n*(axesWidth-2*hMargin);
dy = s(2)*(axesHeight-2*vMargin);

% Start point
switch(lower(location)),
	case 'topleft',
		x0 = axesLeft+hMargin;
		y0 = axesTop-dy-vMargin;
	case 'top',
		x0 = axesLeft+hMargin+dx*n/2;
		y0 = axesTop-dy-vMargin;
	case 'topright',
		x0 = axesRight-hMargin-dx*n;
		y0 = axesTop-dy-vMargin;
	case 'left',
		x0 = axesLeft+hMargin;
		y0 = axesBottom+vMargin+dy/2;
	case 'center',
		x0 = axesLeft+hMargin+dx*n/2;
		y0 = axesBottom+vMargin+dy/2;
	case 'right',
		x0 = axesRight-hMargin-dx*n;
		y0 = axesBottom+vMargin+dy/2;
	case 'bottomleft',
		x0 = axesLeft+hMargin;
		y0 = axesBottom+vMargin;
	case 'bottom',
		x0 = axesLeft+hMargin+dx*n/2;
		y0 = axesBottom+vMargin;
	case 'bottomright',
		x0 = axesRight-hMargin-dx*n;
		y0 = axesBottom+vMargin;
end

% Create all axes
for i = 1:n,
	varargout{i} = axes('position',[x0+(i-1)*dx y0 dx dy]);
end

% Matlab cannot return only varargout, it needs at least one explicit parameter...
a = varargout{1};
varargout = {varargout{2:end}};
