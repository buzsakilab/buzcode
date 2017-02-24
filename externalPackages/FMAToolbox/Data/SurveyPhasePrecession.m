function SurveyPhasePrecession(group,cluster,channel,sessions,varargin)

%SurveyPhasePrecession - Plot phase precession plots for all conditions.
%
%  USAGE
%
%    SurveyPhasePrecession(group,cluster,channel,sessions,<options>)
%
%    group          electrode group (e.g. tetrode ID)
%    cluster        cluster ID
%    channel        channel ID for LFP
%    sessions       optional list of session numbers
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'linearize'   linearization function name or handle
%     'minv'        minimum instantaneous velocity (default = 0)
%     'pixel'       size of the video pixel in cm (no default value)
%     'track'       'linear' or 'circular' (default = 'linear')
%     'figures'     group all plots on one single figure ('single') or create
%                   one figure per session ('multiple', default)
%    =========================================================================
%
%  SEE
%
%    See also PhasePrecession, PlotPhasePrecession.

% Copyright (C) 2009-2012 by Anne Cei, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
	error('Incorrect number of parameters (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).');
end

% Default values
if nargin < 4 || ~isivector(sessions),
	varargin = {sessions,varargin{:}};
	sessions = [];
end
figs = 'multiple';
minv = 0;
track = 'linear';
linearize = GetCustomDefaults('linearize','');

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
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SurveyPhasePrecession">SurveyPhasePrecession</a>'' for details).']);
	end
end

% Do we need to compute phases or can we restore them from a previous run?
phase = Recall('phase',channel);
if isempty(phase);
 lfp = GetLFP(channel);
 fil = FilterLFP(lfp,'passband','theta');
 phase = Phase(fil);
 Store(phase,'phase',channel);
end

% Get and linearize positions, get spikes
positions = GetPositions('coordinates','video','pixel',pixel);
if ~isempty(minv),
	if isempty(pixel),
		error(['Missing pixel size for minimum velocity (type ''help <a href="matlab:help SurveyFiringMaps">SurveyPhasePrecession</a>'' for details).']);
	end
	x = GetPositions('coordinates','real','pixel',pixel);
	v = LinearVelocity(x,30);
	[unused,in] = Threshold(v,'>',minv,'min',10,'max',1);
	positions = positions(in,:);
end
if isempty(linearize),
	warning('No linearization function provided - using position abscissae (ignoring ordinates).');
	x = positions;
else
	x = feval(linearize,positions);
end

% Get spikes and events
spikes = GetSpikes([group cluster]);
start = GetEvents('beginning of .*');
stop = GetEvents('end of .*');
titles = GetEvents('beginning of .*','output','descriptions');

% Sessions
if isempty(sessions),
	sessions = 1:length(start);
end
nSessions = length(sessions);
% Figure information
n = floor(sqrt(nSessions));
m = ceil(nSessions/n);
if strcmp(figs,'single'),
	titleWidth = floor(130/nSessions);
	figureList = figure;
else
	titleWidth = 130;
	figureList = [];
end

% Loop through sessions
Hide('on');
for i = 1:nSessions,
	if strcmp(figs,'single'),
		p = Subpanel(n,m,i);
	else
		figureList = [figureList figure];
		p = Subpanel(1,1,1);
	end
	session = sessions(i);
	% Title
	name = [titles{session}(14:end) ' (' int2str(group) '-' int2str(cluster) ')'];
	SplitTitle(p,name,titleWidth);
	% Restrict data to this session
	sessionSpikes = Restrict(spikes,[start(session) stop(session)]);
	if size(sessionSpikes,1) < 50,
		continue;
	end
	sessionX = Restrict(x,[start(session) stop(session)]);
	if isempty(sessionX), continue; end
	sessionPhase = Restrict(phase,[start(session) stop(session)]);
	[data,stats] = PhasePrecession(sessionX,sessionSpikes,sessionPhase);
	PlotPhasePrecession(data,'parent',p,'track',track);
end
Hide('off');
Hide(figureList,'off');
