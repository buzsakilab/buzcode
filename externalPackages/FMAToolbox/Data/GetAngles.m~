function angles = GetAngles(varargin)

%GetAngles - Get angles (from position samples).
%
%  Angles are returned in radians (in [-pi,pi]). This function requires two
%  or more LEDs (only the first two are used).
%
%  USAGE
%
%    angles = GetAngles(<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'         'all' gets all position samples, 'clean' discards bad
%                    position samples (lights out of boundaries or too close
%                    - see SETTINGS) (default 'clean')
%    =========================================================================

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA SETTINGS;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Default values
mode = 'clean';

if mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help GetAngles">GetAngles</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help GetAngles">GetAngles</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'mode',
			mode = lower(varargin{i+1});
			if ~isstring(mode,'clean','all'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help GetAngles">GetAngles</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GetAngles">GetAngles</a>'' for details).']);
	end
end

angles = [];
positions = DATA.positions;
if isempty(positions), return; end

% Make sure we have at least two LEDs
if size(positions,2) < 5,
	error('Insufficient number of LEDs (type ''help <a href="matlab:help GetAngles">GetAngles</a>'' for details).');
end

% Discard positions with incorrect distance between lights
if strcmp(mode,'clean'),
	distance = sqrt((positions(:,4)-positions(:,2)).^2+(positions(:,5)-positions(:,3)).^2);
	selected = (distance >= SETTINGS.minDistance & distance <= SETTINGS.maxDistance);
	positions = positions(selected,:);
end

% Discard partial/complete undetects (by convention, values of -1)
discard = any(positions(:,2:end) == -1,2);
positions(discard,:) = [];

% Clip if necessary (but issue warning)
m = min(positions(:,2:end));
M = max(positions(:,2:end));
if any(m<0)||any(M(1:2:end)>DATA.maxX|M(2:2:end)>DATA.maxY),
	warning(['Position data should stay within [0,' num2str(DATA.maxX) ']x[0,' num2str(DATA.maxY) ']. The data will now be clipped accordingly.']);
	positions(:,2) = Clip(positions(:,2),0,DATA.maxX);
	positions(:,3) = Clip(positions(:,3),0,DATA.maxY);
end

positions(:,2:end) = Smooth(positions(:,2:end),[5 0]);
angles(:,1) = positions(:,1);
angles(:,2) = angle( positions(:,4)-positions(:,2)+j*(positions(:,5)-positions(:,3)) );
