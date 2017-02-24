function SurveyFiringMaps(group,cluster,varargin)

%SurveyFiringMaps - Plot firing maps for all conditions in the same figure.
%
%  USAGE
%
%    SurveyFiringMaps(group,cluster,<options>)
%
%    group          optional electrode group (e.g. tetrode ID)
%    cluster        optional cluster ID
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'minv'        minimum instantaneous velocity (default = 0)
%     'pixel'       size of the video pixel in cm (no default value)
%    =========================================================================
%
%  SEE
%
%    See also FiringMap, PlotFiringMap.

% Copyright (C) 2009-2011 by Anne Cei, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
units = GetUnits;
minv = [];
pixel = [];

% Which parameters were passed?
if nargin >= 2 && isiscalar(group) && isiscalar(cluster),
	% Group and cluster
	units = units(units(:,1)==group&units(:,2)==cluster,:);
elseif nargin >= 1 && isiscalar(group),
	% Group but not cluster
	units = units(units(:,1)==group,:);
	if nargin >= 2, varargin = {cluster,varargin{:}}; end
else
	% Neither group nor cluster
	switch nargin,
		case 0,
			varargin = {};
		case 1,
			varargin = {group};
		case 2,
			varargin = {group,cluster,varargin{:}};
	end
end
nUnits = size(units,1);

if mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'minv',
			minv = varargin{i+1};
			if ~isdscalar(minv,'>=0'),
				error('Incorrect value for property ''minv'' (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).');
			end
		case 'pixel',
			pixel = varargin{i+1};
			if ~isdscalar(pixel,'>0'),
				error('Incorrect value for property ''pixel'' (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).']);
	end
end

% Get positions
positions = GetPositions;
if isempty(positions), warning('No positions found for current session'); return; end
if ~isempty(minv),
	if isempty(pixel),
		error(['Missing pixel size for minimum velocity (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).']);
	end
	x = GetPositions('coordinates','real','pixel',pixel);
	v = LinearVelocity(x,30);
	[unused,in] = Threshold(v,'>',minv,'min',10,'max',10);
	positions = positions(in,:);
end

% Get start/stop events for each condition
start = GetEvents('beginning of .*');
stop = GetEvents('end of .*');
nSessions = length(start);

titles = GetEvents('beginning of .*','output','descriptions');

% Loop through units
Hide('on');
figureList = [];
for j = 1:nUnits,
	figureList = [figureList figure];
	spikes = GetSpikes(units(j,:));
	% Loop through conditions
	for i = 1:nSessions,
		SquareSubplot(nSessions,i);
		p = Restrict(positions,[start(i) stop(i)]);
		s = Restrict(spikes,[start(i) stop(i)]);
		map = FiringMap(p(:,1:3),s,'nbins',[250 250],'smooth',3,'mintime',0);
		PlotColorMap(map.rate,map.time,'ydir','reverse','bar','off','cutoffs',[0 10]);
		title = titles{i}(14:end);
		SplitTitle([title ' (' int2str(units(j,1)) '-' int2str(units(j,2)) ')'],round(150/sqrt(nSessions)));
	end
end
Hide('off');
Hide(figureList,'off');