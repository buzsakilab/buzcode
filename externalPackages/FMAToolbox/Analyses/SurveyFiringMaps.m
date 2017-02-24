function [maps,stats] = SurveyFiringMaps(units,varargin)

%SurveyFiringMaps - Compute and plot firing maps for all subsessions.
%
%  USAGE
%
%    [maps,stats] = SurveyFiringMaps(units,<options>)
%
%    units          optional list of units, i.e. [electrode group, cluster]
%                   pairs; set cluster to -1 to process all clusters
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'minv'        minimum instantaneous velocity (default = 0)
%     'pixel'       size of the video pixel in cm (no default value)
%     'show'        set to 'off' to compute but not plot data (default = 'on')
%    =========================================================================
%
%  OUTPUT
%
%    The outputs are the same as for <a href="matlab:help MapStats">MapStats</a>, except for map.z which is
%    replaced by map.rate.
%
%  SEE
%
%    See also FiringMap, PlotColorMap.

% Copyright (C) 2009-2013 by Anne Cei, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
minv = [];
pixel = [];
show = 'on';

% No unit list provided?
if ~(nargin >= 1 && isimatrix(units)),
	units = GetUnits;
	varargin = {units,varargin{:}};
else
	all = units(:,2) == -1;
	groups = unique(units(all,1));
	units(all,:) = [];
	units = unique([units;GetUnits(groups)],'rows');
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
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).']);
	end
end

% Get positions
positions = GetPositions;
positions=[positions(:,1) positions(:,4:5)];
positions(234:244,:)
if isempty(positions), warning('No positions found for current subsession'); return; end
if ~isempty(minv),
	if isempty(pixel),
		error(['Missing pixel size for minimum velocity (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).']);
	end
	x = GetPositions('coordinates','real','pixel',pixel);
	v = LinearVelocity(x,30);
	[~,in] = Threshold(v,'>',minv,'min',10,'max',2);
	positions = positions(in,:);
end

% Get start/stop events for each subsession
start = GetEvents('beginning of .*');
stop = GetEvents('end of .*');
nSubsessions = length(start);

if strcmp(show,'on'),
	titles = GetEventTypes('beginning of .*');
	status = Hide('status');
	Hide('on');
	figureList = [];
end

try
	n = 0;
	% Loop through units
	for j = 1:nUnits,
		if strcmp(show,'on'), figureList = [figureList figure]; end
		spikes = GetSpikes(units(j,:));
		% Loop through subsessions
		for i = 1:nSubsessions,
			p = Restrict(positions,[start(i) stop(i)]);
			s = Restrict(spikes,[start(i) stop(i)]);
			if nargin == 0,
				map = FiringMap(p(:,1:3),s,'nbins',[250 250],'smooth',5,'mintime',0);
			else
				[map,st] = FiringMap(p(:,1:3),s,'nbins',[250 250],'smooth',5,'mintime',0);
				st.unit = units(j,:);
				st.subsession = i;
				n = n + 1;
				stats(n) = st;
				maps(n) = map;
			end
			if strcmp(show,'on'),
				SquareSubplot(nSubsessions,i);
				PlotColorMap(map.rate,map.time,'ydir','reverse','bar','off','cutoffs',[0 10],'gamma',2);
				title = titles{i}(14:end);
				SplitTitle([title ' (' int2str(units(j,1)) '-' int2str(units(j,2)) ')'],round(150/sqrt(nSubsessions)));
			end
		end
	end
catch err
	if strcmp(show,'on'),
		Hide(status);
		Hide(figureList,status);
	end
	err.rethrow;
end

if strcmp(show,'on'),
	Hide(status);
	Hide(figureList,status);
end