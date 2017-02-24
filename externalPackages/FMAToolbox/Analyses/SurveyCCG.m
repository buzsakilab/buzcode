function [data,t] = SurveyCCG(units,varargin)

%SurveyCCG - Compute and plot cross-correlograms (CCGs) for all subsessions.
%
%  USAGE
%
%    [data,t] = SurveyCCG(units,<options>)
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
%     'duration'    duration in s of each xcorrelogram (default = 1)
%    =========================================================================
%
%  OUTPUT
%
%    The outputs are the same as for <a href="matlab:help CCG">CCG</a>.
%
%  SEE
%
%    See also CCG, SurveyFiringMaps, SurveyPhasePrecession.

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
duration = 1;

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
	error('Incorrect number of parameters (type ''help <a href="matlab:help SurveyCCG">SurveyCCG</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help SurveyCCG">SurveyCCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'minv',
			minv = varargin{i+1};
			if ~isdscalar(minv,'>=0'),
				error('Incorrect value for property ''minv'' (type ''help <a href="matlab:help SurveyCCG">SurveyCCG</a>'' for details).');
			end
		case 'pixel',
			pixel = varargin{i+1};
			if ~isdscalar(pixel,'>0'),
				error('Incorrect value for property ''pixel'' (type ''help <a href="matlab:help SurveyCCG">SurveyCCG</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help SurveyCCG">SurveyCCG</a>'' for details).');
			end
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help SurveyCCG">SurveyCCG</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SurveyCCG">SurveyCCG</a>'' for details).']);
	end
end

% Get spikes (with unit IDs)
spikes = GetSpikes(units,'output','full');

% Get positions, compute velocity and discard spikes during slow periods
positions = GetPositions;
if isempty(positions), warning('No positions found for current session'); return; end
if ~isempty(minv),
	if isempty(pixel),
		error(['Missing pixel size for minimum velocity (type ''help <a href="matlab:help SurveyCCG">SurveyCCG</a>'' for details).']);
	end
	x = GetPositions('coordinates','real','pixel',pixel);
	x = Interpolate(x,spikes(:,1));
	v = LinearVelocity(x,30);
	[~,in] = Threshold(v,'>',minv,'min',10,'max',10);
	spikes = spikes(in,:);
end

% Get start/stop events for each subsession
start = GetEvents('beginning of .*');
stop = GetEvents('end of .*');
nSubsessions = length(start);

% Number units from 1 to n
n = size(units,1);
groups = ones(size(spikes(:,1)));
for i = 1:n,
	this = spikes(:,2)==units(i,1)&spikes(:,3)==units(i,2);
	groups(this) = i;
end
spikes = spikes(:,1);

if strcmp(show,'on'),
	titles = GetEventTypes('beginning of .*');
	status = Hide('status');
	Hide('on');
	figureList = [];
end

try
	% Loop through subsessions
	for i = 1:nSubsessions,
		if strcmp(show,'on'), figureList = [figureList figure]; end
		in = InIntervals(spikes,[start(i) stop(i)]);
		s = spikes(in);
		g = groups(in);
		[ccg,t] = CCG(s,g,'binSize',duration/30,'duration',duration);
		d.ccg = ccg;
		d.subsession = i;
		data(i) = d;
		if strcmp(show,'on'),
			PlotCCG(t,ccg);
			subplot(n,n,n*n);
			axis off;
			title(titles{i}(14:end));
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