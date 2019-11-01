function positions = GetPositions(varargin)

%GetPositions - Get position samples.
%
%  USAGE
%
%    positions = GetPositions(<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'         'all' gets all position samples, 'clean' discards bad
%                    position samples (lights undetected, out of boundaries,
%                    too close, etc. - see SETTINGS) (default 'clean')
%     'coordinates'  'video' gets position samples in video frame coordinates
%                    (e.g., [0..368]x[0..284]), 'normalized' uses normalized
%                    values (i.e. [0..1]x[0..1]), 'real' converts pixels to
%                    centimeters (default 'normalized')
%     'pixel'        size of the video pixel in cm (no default value)
%     'discard'      discard position samples when one or more LED is missing
%                    ('partial', default) or only when no LEDs are detected
%                    ('none', coordinates for missing LEDs are set to NaN)
%     'distances'    min and max distances [m M] allowed between LEDs
%                    (in pixels, default [0 Inf])
%    =========================================================================
%
%  EXAMPLES
%
%    p = GetPositions;
%    p = GetPositions('mode','all');
%
%  CUSTOM DEFAULTS
%
%    Properties 'mode', 'coordinates', 'pixel' and 'distances' can have custom
%    default values (type 'help <a href="matlab:help CustomDefaults">CustomDefaults</a>' for details).

% Copyright (C) 2004-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA SETTINGS;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Default values (customizable defaults must be empty at this point)
mode = '';
coordinates = '';
pixel = [];
distances = [];
discard = 'partial';

if mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'mode',
			mode = lower(varargin{i+1});
			if ~isstring_FMAT(mode,'clean','all'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).');
			end
		case 'coordinates',
			coordinates = lower(varargin{i+1});
			if ~isstring_FMAT(coordinates,'video','normalized','real'),
				error('Incorrect value for property ''coordinates'' (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).');
			end
		case 'pixel',
			pixel = varargin{i+1};
			if ~isdscalar(pixel,'>0'),
				error('Incorrect value for property ''pixel'' (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).');
			end
		case 'discard',
			discard = lower(varargin{i+1});
			if ~isstring_FMAT(discard,'partial','none'),
				error('Incorrect value for property ''discard'' (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).');
			end
		case 'distances',
			distances = lower(varargin{i+1});
			if ~isdvector(distances,'#2','>=0'),
				error('Incorrect value for property ''distances'' (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).']);
	end
end

% Customizable defaults (if not already set)
mode = GetCustomDefaults(mode,'mode','clean');
coordinates = GetCustomDefaults(coordinates,'coordinates','normalized');
pixel = GetCustomDefaults(pixel,'pixel',[]);
distances = GetCustomDefaults(distances,'distances',[0 Inf]);

minDistance = distances(1);
maxDistance = distances(2);

if strcmp(coordinates,'real') & isempty(pixel),
	error(['Missing ''pixel'' property-value pair (type ''help <a href="matlab:help GetPositions">GetPositions</a>'' for details).']);
end

positions = DATA.positions;
if isempty(positions), return; end

% Discard positions with incorrect distance between lights
if strcmp(mode,'clean'),
	if size(positions,2) == 5,
		% Two head lamps
		distance = sqrt((positions(:,4)-positions(:,2)).^2+(positions(:,5)-positions(:,3)).^2);
		selected = ( ...
				distance >= minDistance ...
				& distance <= maxDistance ...
			);
		positions = positions(selected,:);
	end
end

% Mark undetects (by convention, values of -1)
if strcmp(discard,'none'),
	undetected = all(positions(:,2:end) == -1,2);
	partiallyDetected = positions == -1;
	positions(partiallyDetected) = NaN;
else
    undetected = any(positions(:,2:end) == -1,2);

end
%  positions(undetected,2:end) = -1; % Make sure all values are -1
good = ~undetected;

% Clip if necessary (but issue warning)
m = min(positions(good,2:end));
M = max(positions(good,2:end));
%  if any(m<0)||any(M(1:2:end)>DATA.maxX|M(2:2:end)>DATA.maxY),
%  	warning(['Position data should stay within [0,' num2str(DATA.maxX) ']x[0,' num2str(DATA.maxY) ']. The data will now be clipped accordingly.']);
%  	positions(good,2) = Clip(positions(good,2),0,DATA.maxX);
%  	positions(good,3) = Clip(positions(good,3),0,DATA.maxY);
%  end

if strcmp(coordinates,'normalized'),
	% Normalize coordinates
	maxima = M(1:2);%[DATA.maxX DATA.maxY];
	positions(good,2:end) = positions(good,2:end) ./ repmat(maxima,sum(good),(size(positions,2)-1)/2);
elseif strcmp(coordinates,'real'),
	% Convert to cm
	positions(good,2:end) = positions(good,2:end) * pixel;
end

% Discard undetects
if strcmp(mode,'clean'),
	positions(undetected,:) = [];
end
