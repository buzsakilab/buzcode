function m = Bright(n,varargin)

%Bright - Bright colormap (similar to HSV or JET, but brighter).
%
%  USAGE
%
%    m = Bright(n,<options>)
%
%    n              optional number of rows in output matrix (default = 100)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'hgamma'      gamma-like correction for hue (1 = no correction, default)
%     'type'        either 'linear' or 'circular' (default 'linear')
%     'stops'       hue stops (default [2/3 0] linear, [0 1] circular)
%    =========================================================================

% Copyright (C) 2009-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
hgamma = 1;
type = 'linear';
stops = [];
linearStops = [2/3 0];
circularStops = [0 1];

% Optional parameter
if nargin < 1,
	n = 100;
elseif ischar(n),
	varargin = {n,varargin{:}};
	n = 100;
elseif ~isdscalar(n,'>=0'),
	error('Incorrect value for ''n'' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Bright">Bright</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'hgamma',
			hgamma = varargin{i+1};
			if ~isdscalar(hgamma,'>=0'),
			error('Incorrect value for property ''hgamma'' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).');
			end
		case 'stops',
			stops = varargin{i+1};
			if ~isdvector(stops,'>=0','<=1') || length(stops) < 2,
				error('Incorrect value for property ''stops'' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).');
			end
		case 'type',
			type = lower(varargin{i+1});
			if ~isstring_FMAT(type,'linear','circular'),
			error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).']);
	end
end

% Default stops
if isempty(stops),
	if strcmp(type,'circular'),
		stops = circularStops;
	else
		stops = linearStops;
	end
end
if strcmp(type,'circular') && stops(1) ~= stops(end) && abs(stops(1)-stops(end)) ~= 1,
	stops(end+1) = stops(1);
end

% Number of color band (transitions between stop pairs)
nBands = length(stops)-1;

% Construct color bands
hsv = [];
nn = round(n/nBands);
for i = 1:nBands,
	hsv = [hsv ; linspace(0,1,nn)'.^(1/hgamma)*(stops(i+1)-stops(i))+stops(i)];
end
n = length(hsv);

% Set saturation and value, then convert to RGB
hsv(:,2) = 1;
hsv(:,3) = 1;
m = hsv2rgb(hsv);

