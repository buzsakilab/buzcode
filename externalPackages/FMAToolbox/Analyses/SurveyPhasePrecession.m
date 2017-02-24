function [data,stats] = SurveyPhasePrecession(group,clusters,channel,subsessions,varargin)

%SurveyPhasePrecession - Compute and plot phase precession plots for all subsessions.
%
%  USAGE
%
%    [data,stats] = SurveyPhasePrecession(group,clusters,channel,subsessions,<options>)
%
%    group          electrode group (e.g. tetrode ID)
%    clusters       list of cluster IDs
%    channel        channel ID for LFP
%    subsessions    optional list of subsession numbers
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'linearize'   linearization function name or handle
%     'minv'        minimum instantaneous velocity (default = 0)
%     'pixel'       size of the video pixel in cm (no default value)
%     'track'       'linear' or 'circular' (default = 'linear')
%     'show'        set to 'off' to compute but not plot data (default = 'on')
%     'figures'     for each cluster, group all plots on one single figure
%                   ('single') or create one figure per subsession ('multiple',
%                   default)
%     'slopes'      search start value for slopes (default = 0) for each
%                   subsession
%    =========================================================================
%
%  OUTPUT
%
%    The outputs are the same as for <a href="matlab:help PhasePrecession">PhasePrecession</a>.
%
%  CUSTOM DEFAULTS
%
%    Poperties 'linearize', 'pixel' and 'minv' can have custom default values
%    (type 'help <a href="matlab:help CustomDefaults">CustomDefaults</a>' for details).
%
%  SEE
%
%    See also PhasePrecession, PlotPhasePrecession.

% Copyright (C) 2009-2013 by Anne Cei, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
	error('Incorrect number of parameters (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');
end

% Default values
if nargin < 4 || ~isivector(subsessions),
	if nargin >= 4, varargin = {subsessions,varargin{:}}; end
	subsessions = [];
end
figs = 'multiple';
minv = 0;
track = 'linear';
linearize = '';
show = 'on';
slopes = [];
nBins = 200;

% Unit list
if clusters == -1,
	units = GetUnits(group);
	clusters = units(:,2);
end

if mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'linearize',
			linearize = lower(varargin{i+1});
			if isa(linearize,'function_handle'),
				linearize = func2str(linearize);
			end
			if ~ischar(linearize) || isempty(which(linearize)),
				error('Incorrect value for property ''linearize'' (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');
			end
		case 'minv',
			minv = varargin{i+1};
			if ~isdscalar(minv,'>=0'),
				error('Incorrect value for property ''minv'' (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');
			end
		case 'pixel',
			pixel = varargin{i+1};
			if ~isdscalar(pixel,'>0'),
				error('Incorrect value for property ''pixel'' (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');
			end
		case 'figures',
			figs = lower(varargin{i+1});
			if ~isstring_FMAT(figs,'single','multiple'),
				error('Incorrect value for property ''figs'' (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');
			end
		case 'track',
			track = lower(varargin{i+1});
			if ~isstring_FMAT(track,'linear','circular'),
				error('Incorrect value for property ''track'' (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');		
			end
		case 'slopes',
			slopes = varargin{i+1};
			if ~isdvector(slopes),
				error('Incorrect value for property ''slopes'' (type ''help <a href="matlab:help CircularRegression">PhasePrecession</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).']);
	end
end

% Default linearization function
if isempty(linearize),
	linearize = GetCustomDefaults('linearize','');
end
% Default pixel size
if isempty(pixel),
	linearize = GetCustomDefaults('pixel','');
end
% Default minimum velocity
if isempty(minv),
	linearize = GetCustomDefaults('minv','');
end


% Do we need to compute phases or can we restore them from a previous run?
phase = Recall('phase',channel);
if isempty(phase);
	lfp = GetLFP(channel);
	fil = FilterLFP(lfp,'passband','theta');
	phase = Phase(fil);
	Store(phase,'phase',channel);
end

% Get and linearize positions,
positions = GetPositions;
if isempty(linearize),
	warning('No linearization function provided - using position abscissae (ignoring ordinates).');
	x = positions(:,1:2);
else
	x = feval(linearize,positions);
end
if minv ~= 0,
	if isempty(pixel),
		error(['Missing pixel size for minimum velocity (type ''help <a href="matlab:help SurveyFiringMaps">SurveyPhasePrecession</a>'' for details).']);
	end
	v = LinearVelocity(GetPositions('coordinates','real','pixel',pixel),30);
	[~,in] = Threshold(v,'>',minv,'min',10,'max',1);
	x = x(in,:);
end

% Get events
start = GetEvents('beginning of .*');
stop = GetEvents('end of .*');

% Sessions
if isempty(subsessions),
	subsessions = 1:length(start);
end
nSubsessions = length(subsessions);
if isempty(slopes),
	slopes = zeros(size(subsessions));
elseif length(slopes) ~= length(subsessions),
	error('Slopes and subsessions do not have the same length (type ''help <a href="matlab:help CircularRegression">PhasePrecession</a>'' for details).');
end

% Figure information
if strcmp(show,'on'),
	titles = GetEventTypes('beginning of .*');
	n = floor(sqrt(nSubsessions));
	m = ceil(nSubsessions/n);
	if strcmp(figs,'single'),
		titleWidth = floor(130/nSubsessions);
	else
		titleWidth = 130;
	end
	figureList = [];
	status = Hide('status');
	Hide('on');
end

try
	k = 0;
	% Loop through units
	for c = 1:length(clusters),
		% Get spikes for this cluster
		cluster = clusters(c);
		spikes = GetSpikes([group cluster]);
		if strcmp(show,'on') && strcmp(figs,'single'),
			figureList = [figureList figure];
		end
		
		% Loop through subsessions
		for i = 1:nSubsessions,
			subsession = subsessions(i);
			% Restrict data to this subsession
			subsessionSpikes = Restrict(spikes,[start(subsession) stop(subsession)]);
			if size(subsessionSpikes,1) < 50,
				continue;
			end
			subsessionX = Restrict(x,[start(subsession) stop(subsession)]);
			if isempty(subsessionX), continue; end
			
	%  f=gcf;figureList = [figureList figure];PlotXY(subsessionX);figure(f);
			subsessionPhase = Restrict(phase,[start(subsession) stop(subsession)]);
			[d,st] = PhasePrecession(subsessionX,subsessionSpikes,subsessionPhase,'slope',slopes(i));
			st.unit = [group cluster];
			st.subsession = i;
			st.channel = channel;
			k = k + 1;
			data(k) = d;
			stats(k) = st;
			if strcmp(show,'on'),
				if strcmp(figs,'single'),
					p = Subpanel(n,m,i);
				else
					figureList = [figureList figure];
					p = Subpanel(1,1,1);
				end
				% Title
				name = [titles{subsession}(14:end) ' (' int2str(group) '-' int2str(cluster) ')'];
				SplitTitle(p,name,titleWidth);
				PlotPhasePrecession(d,st,'parent',p,'track',track);
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