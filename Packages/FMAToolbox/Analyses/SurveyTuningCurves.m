function [curves,stats] = SurveyTuningCurves(units,varargin)

%SurveyTuningCurves - Compute and plot tuning curves for all subsessions.
%
%  USAGE
%
%    [curves,stats] = SurveyTuningCurves(units,<options>)
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
%     'type'        'cartesian' (default) or 'polar' plot
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
type = 'cartesian';

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
		case 'type',
			type = varargin{i+1};
			if ~isstring_FMAT(type,'polar','cartesian'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).']);
	end
end

% Get positions
positions = GetPositions;
if isempty(positions), warning('No positions found for current subsession'); return; end
if ~isempty(minv),
	if isempty(pixel),
		error('Missing pixel size for minimum velocity (type ''help <a href="matlab:help SurveyFiringMaps">SurveyFiringMaps</a>'' for details).');
	end
	x = GetPositions('coordinates','real','pixel',pixel);
	v = LinearVelocity(x,30);
	[~,in] = Threshold(v,'>',minv,'min',10,'max',2);
	positions = positions(in,:);

end

% Get angles and normalize in [0,1]
angles = GetAngles;
angles(:,2) = (angles(:,2)+pi)/(2*pi);

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
			a = Restrict(angles,[start(i) stop(i)]);
			s = Restrict(spikes,[start(i) stop(i)]);
			if nargout == 0,
				curve = FiringCurve(a,s,'nbins',120,'mintime',0,'type','cl');
			else
				[curve,st] = FiringCurve(a,s,'nbins',120,'mintime',0,'type','cl');
				st.unit = units(j,:);
				st.subsession = i;
				n = n + 1;
				stats(n) = st;
				curves(n) = curve;
			end
			if strcmp(show,'on'),
				SquareSubplot(nSubsessions,i);
				if strcmp(type,'cartesian'),
					plot(curve.x,curve.rate);
				else
					polar(curve.x*2*pi,curve.rate);
				end
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