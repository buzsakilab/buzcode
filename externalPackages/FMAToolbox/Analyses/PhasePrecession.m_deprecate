function [data,stats] = PhasePrecession(positions,spikes,phases,varargin)

%PhasePrecession - Compute spike phase precession.
%
% Compute spike phase precession using the methods of O'Keefe and Recce (1993),
% i.e. spike phase vs position, and Harris et al. (2002), i.e. spike phase vs
% spike rate. Single-lap phase precession is also computed.
%
%  USAGE
%
%    [data,stats] = PhasePrecession(positions,spikes,phases,<options>)
%
%    positions      linearized position samples (normalized to [0,1])
%    spikes         spike timestamps
%    phases         phase samples (see <a href="matlab:help Phase">Phase</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'maxGap'      time gaps between successive position samples exceeding
%                   this threshold (e.g. undetects) will not be interpolated
%                   (default = 100 ms)
%     'boundaries'  onset and offset for single-lap phase precession can be
%                   determined either automatically based on spike count
%                   ('count', default) or using explicit firing field
%                   boundaries ([Xstart Xstop], where each X is in [0..1])
%    =========================================================================
%
%  NOTE
%
%    To compute only phase vs rate, positions are not required and can be left
%    empty, e.g. PhasePrecession([],spikes,phases).
%
%  OUTPUT
%
%    data.x           position samples
%    data.position    spike phase vs position:
%     .x               x coordinates at each spike timestamp
%     .t               spike times
%     .phase           spike phases (in radians)
%     .lap             lap number
%
%    data.rate        spike phase vs spike rate:
%     .r               spike rates at each spike timestamp
%     .t               spike times
%     .phase           spike phases (in radians)
%     .lap             lap number
%
%    Additional statistics computed using phase vs position data:
%
%    stats.slope      phase precession slope (via circular regression)
%    stats.intercept  phase precession intercept (via circular regression)
%    stats.r2         coefficient of determination for circular regression
%    stats.center     center x for early, middle and late subfields
%    stats.mean       mean phase for early, middle and late subfields
%    stats.var        phase variance for early, middle and late subfields
%    stats.std        phase standard deviation for early, middle and late subfields
%    stats.conf       95% confidence intervals
%    stats.all        cell array containing all phases for each subfield
%                     (useful for population analyses)
%
%    stats.lap        single-lap statistics
%     .slope          phase precession slope (via circular regression)
%     .intercept      phase precession intercept (via circular regression)
%     .r2             coefficient of determination for circular regression
%
%  SEE
%
%    See also Phase, PlotPhasePrecession.

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
maxGap = 0.1;
boundaries = 'count';

% Check number of parameters
if nargin < 3 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
end

% Check parameter sizes
if ~isempty(positions) && size(positions,2) < 2,
	error('Parameter ''positions'' should have at least 2 columns (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
end
if ~isdvector(spikes),
	error('Parameter ''spikes'' is not a vector (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
end
if size(phases,2) ~= 2,
	error('Parameter ''phases'' is not a Nx2 matrix (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
end
isradians(phases(:,2));

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'boundaries',
			boundaries = lower(varargin{i+1});
			if ~isdvector(boundaries,'#2','>=0','<=1') && ~isstring_FMAT(boundaries,'count'),
				error('Incorrect value for property ''boundaries'' (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
			end
		case 'maxgap',
			maxGap = varargin{i+1};
			if ~isdscalar(maxGap,'>0'),
				error('Incorrect value for property ''maxGap'' (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).']);
	end
end

% Default values
data.x = positions;
data.position.x = [];
data.position.t = [];
data.position.phase = [];
data.position.lap = [];
data.rate.r = [];
data.rate.t = [];
data.rate.phase = [];
data.rate.lap = [];
stats.slope = [];
stats.r2 = [];
stats.center = [0 0 0];
stats.mean = [0 0 0];
stats.var = [0 0 0];
stats.std = [0 0 0];
stats.conf = [0 0 0;0 0 0];

% Compute spike phases
if isempty(spikes), return; end
spikePhases = Interpolate(phases,spikes,'trim','off','type','circular');
if isempty(spikePhases), return; end

% Interpolate positions at spike times
if ~isempty(positions),
	% Make sure positions are normalized
	if max(positions(:,2)) > 1 || min(positions(:,2)) < 0,
		positions(:,2) = ZeroToOne(positions(:,2));
		warning('Parameter ''positions'' should contain values in [0 1]. The data will now be transformed accordingly.');
	end
	[x,ignored] = Interpolate(positions,spikes,'trim','off','maxGap',maxGap);
	data.position.t = spikes(~ignored);
	data.position.x = x(:,2);
	data.position.phase = spikePhases(~ignored,2);
end

% Compute unwrapped spike phases
unwrapped = [phases(:,1) unwrap(phases(:,2))];
[spikeUnwrappedPhases,ignored] = Interpolate(unwrapped,spikes,'trim','off');
% Get beginning and end unwrapped phases of two cycles surrounding each spike
startUnwrappedPhases = spikeUnwrappedPhases(:,2)-2*pi;
stopUnwrappedPhases = spikeUnwrappedPhases(:,2)+2*pi;

%%%%
%  figure;PlotXY(unwrapped);hold on;
%  PlotHVLines(spikes,'v','k');PlotXY(spikeUnwrappedPhases,'r+');PlotHVLines(startUnwrappedPhases,'h','g');PlotHVLines(stopUnwrappedPhases,'h','r');
%  figure;PlotXY(phases);hold on;
%  PlotTicks(spikes,pi,'k');
%%%%

% Compute instantaneous rate at spike times
data.rate.r = CountInIntervals(spikeUnwrappedPhases(:,2),[startUnwrappedPhases stopUnwrappedPhases]);
data.rate.t = spikes;
data.rate.phase = spikePhases(:,2);


% Determine lap number for each spike (for single-lap phase precession)

if isdvector(boundaries,'#2') && ~isempty(positions),

	% Firing curve boundaries
	fieldStart = boundaries(1);
	fieldStop = boundaries(2);
	fieldSize = abs(fieldStart-fieldStop);

	% Let x1 be current position, x0 be previous position and x2 be next position
	x0 = [positions(1,2);positions(1:end-1,2)];
	x1 = positions(:,2);
	x2 = [positions(2:end,2);positions(end,2)];

	% Determine when the animal enters/leaves the field:
	% 1) the animal enters the field when x1 is inside and x0 is not inside and near x1
	% 2) the animal leaves the field when x1 is inside and x2 is not inside and near x1
	if fieldStart < fieldStop,
		x0_inside = x0 >= fieldStart & x0 <= fieldStop;
		x1_inside = x1 >= fieldStart & x1 <= fieldStop;
		x2_inside = x2 >= fieldStart & x2 <= fieldStop;
	else
		x0_inside = x0 >= fieldStart | x0 <= fieldStop;
		x1_inside = x1 >= fieldStart | x1 <= fieldStop;
		x2_inside = x2 >= fieldStart | x2 <= fieldStop;
	end
	x0_near_x1 = abs(x0-x1)<0.1*fieldSize;
	x2_near_x1 = abs(x2-x1)<0.1*fieldSize;
	start = x1_inside & ~(x0_inside & x0_near_x1);
	stop = x1_inside & ~(x2_inside & x2_near_x1);
	start = positions(start,1);
	stop = positions(stop,1);

	if any(start) && any(stop),
		% Get rid of duplicate events (e.g. the rat enters twice before it exits)
		s = [start ones(size(start));stop 2*ones(size(stop))];
		s = sortrows(s);
		ds = [1;diff(s(:,2))];
		s = s(ds~=0,:);
		start = s(s(:,2)==1);
		stop = s(s(:,2)==2);
		% Drop incomplete traversals of the field
		if start(1) > stop(1), stop(1) = []; end
		if start(end) > stop(end), start(end) = []; end
		traversals = [start stop];

		% Determine lap number
		[unused,lap] = InIntervals(data.position.t,traversals);
		data.position.lap = lap;
		[unused,lap] = InIntervals(data.rate.t,traversals);
		data.rate.lap = lap;
	end

elseif strcmp(boundaries,'count'),

	% Get beginning and end unwrapped phases of eight cycles before each spike
	startUnwrappedPhases = spikeUnwrappedPhases(:,2)-8*2*pi;
	stopUnwrappedPhases = spikeUnwrappedPhases(:,2);
	rateBefore = CountInIntervals(spikeUnwrappedPhases,[startUnwrappedPhases stopUnwrappedPhases]);
	% Get beginning and end unwrapped phases of eight cycles after each spike
	startUnwrappedPhases = spikeUnwrappedPhases(:,2);
	stopUnwrappedPhases = spikeUnwrappedPhases(:,2)+8*2*pi;
	rateAfter = CountInIntervals(spikeUnwrappedPhases,[startUnwrappedPhases stopUnwrappedPhases]);
	% Find periods of minimal (resp. maximal) firing rates
	minimalBefore = rateBefore <= 4;
	maximalAfter = rateAfter >= 16;
	maximalBefore = rateBefore >= 16;
	minimalAfter = rateAfter <= 4;
	start = minimalBefore & maximalAfter;
	stop = maximalBefore & minimalAfter;
	% Drop incomplete traversals of the field
	if any(start) && any(stop),
		if start(1) >= stop(1), stop(1) = []; end
		if start(end) >= stop(end), start(end) = []; end
		traversals = [start stop];
		% Determine lap number
		[unused,lap] = InIntervals(data.position.t,traversals);
		data.position.lap = lap;
		[unused,lap] = InIntervals(data.rate.t,traversals);
		data.rate.lap = lap;
	end

end

% Statistics: regression + circular means and variances for early, middle and late subfields
if ~isempty(positions) && nargout >= 2,
	ok = ~isnan(data.position.x);
	[beta,r2] = CircularRegression(data.position.x(ok),data.position.phase(ok));
	stats.slope = beta(1);
	stats.intercept = beta(2);
	stats.r2 = r2;
	for lap = 1:max(data.rate.lap),
		thisLap = data.position.lap == lap;
		[beta,r2] = CircularRegression(data.position.x(thisLap&ok),data.position.phase(thisLap&ok),stats.slope);
		stats.lap.slope(lap,1) = beta(1);
		stats.lap.intercept(lap,1) = beta(2);
		stats.lap.r2(lap,1) = r2;
	end
	x0 = min(data.position.x);
	x1 = max(data.position.x);
	dx = x1-x0;
	for i = 1:3,
		ok = data.position.x > x0+dx/3*(i-1) & data.position.x < x0+dx/3*i;
		if sum(ok) == 0, continue;	end
		stats.all{i} = data.position.phase(ok);
		stats.x(i) = x0+1.5*dx/3*(i-1);
		[stats.var(i),stats.std(i)] = CircularVariance(data.position.phase(ok));
		[stats.mean(i),stats.conf(:,i)] = CircularConfidenceIntervals(data.position.phase(ok));
	end
end

