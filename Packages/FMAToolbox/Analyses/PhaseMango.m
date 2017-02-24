function [map,data] = PhaseMango(spikes,phases,varargin)

%PhaseMango - Compute phase as a function of spike rate and acceleration.
%
% Compute spike phase as a function of spike rate and acceleration, as in
% Harris et al. (2002). The name of the function refers to the aspect of
% the resulting phase plot.
%
%  USAGE
%
%    [map,data] = PhaseMango(spikes,phases,<options>)
%
%    spikes         spike timestamps
%    phases         phase samples (see <a href="matlab:help Phase">Phase</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'       number of horizontal and vertical bins (default = [50 50])
%     'limits'      [maxRate minAccel maxAccel] (default = [8 -2 2])
%     'maxGap'      time gaps between successive position samples exceeding
%                   this threshold (e.g. undetects) will not be interpolated
%                   (default = 100 ms)
%     'boundaries'  onset and offset for single-lap phase precession can be
%                   determined either automatically based on spike count
%                   ('count', default) or using explicit firing field
%                   boundaries ([Xstart Xstop], where each X is in [0..1])
%    =========================================================================
%
%  OUTPUT
%
%    map.x                x bins
%    map.y                y bins
%    map.phase            phase map
%    map.count            count map
%
%    data.rate            spike rate for each cycle
%    data.acceleration    spike acceleration for each cycle
%    data.phase           mean spike phase (in radians) for each cycle
%
%  SEE
%
%    See also Phase, PhasePrecession, PlotPhaseMango.

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
maxGap = 0.1;
smooth = 2;
nBins = 50;
limits = [8 -2 2];

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).');
end

% Check parameter sizes
if ~isdvector(spikes),
	error('Parameter ''spikes'' is not a vector (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).');
end
if size(phases,2) ~= 2,
	error('Parameter ''phases'' is not a Nx2 matrix (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).');
end
isradians(phases(:,2));

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdvector(smooth,'>=0') || length(smooth) > 2,
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).');
			end
		case 'nbins',
			nBins = varargin{i+1};
			if ~isivector(nBins,'>0') || length(nBins) > 2,
				error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).');
			end
		case 'limits',
			limits = varargin{i+1};
			if ~isivector(limits,'#3') || limits(1) <0 || limits(3) < limits(2),
				error('Incorrect value for property ''limits'' (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).');
			end
		case 'maxgap',
			maxGap = varargin{i+1};
			if ~isdscalar(maxGap,'>0'),
				error('Incorrect value for property ''maxGap'' (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PhaseMango">PhaseMango</a>'' for details).']);
	end
end

% Default values
data.time = [];
data.rate = [];
data.acceleration = [];
data.phase = [];
map.x = [];
map.y = [];
map.phase = [];
map.count = [];

% Compute spike phases
if isempty(spikes), return; end
spikePhases = Interpolate(phases,spikes,'trim','off','type','circular');
spikePhases = exp(j*spikePhases(:,2));
if isempty(spikePhases), return; end

% Determine theta cycles (start/stop)
phases(:,2) = wrap(phases(:,2),1); % in [-pi,pi]
start = phases(ZeroCrossings(phases),1);
stop = start(2:end);
start = start(1:end-1);
data.time = start;

% Determine to which theta cycle each spike belongs
[~,spikeCycles] = InIntervals(spikes,[start stop]);
partialCycles = spikeCycles==0;
spikeCycles(partialCycles) = [];

% Compute number of spikes in each cycle
nSpikes = Accumulate(spikeCycles,1,size(start));

% Compute average phase in each cycle
spikePhases(partialCycles) = [];
data.phase = angle(Accumulate(spikeCycles,spikePhases,size(start))./nSpikes);

phase = Accumulate(spikeCycles,spikePhases,size(start))./nSpikes;
phase(nSpikes==0) = 0;
kernel = ones(7,1)/7;
data.phase = angle(conv(phase,kernel,'same'));

% Estimate spike rate and acceleration in each cycle

% Rate = average number of spikes over 7 surrounding cycles, which can be
% efficiently computed as the convolution of nSpikes with the kernel (1/7,..,1/7)
kernel = ones(7,1)/7;
data.rate = conv(nSpikes,kernel,'same');
% Acceleration = slope of linear regression over 7 surrounding cycles
% The general formula is b = (E[xy]-E[x]E[y]) / (E[x²]-E[x]²)
% To compute this efficiently, we first note that the slope is unchanged
% by translations along x the axis: we choose x in {-3..3}, so that E[x]=0.
% We then compute E[x²] = 4, and end up with b = 1/4 E[xy], which can be
% computed as the convolution of y with the kernel (3,2,1,0,-1,-2,-3)/7.
kernel = (3:-1:-3)/7;
data.acceleration = conv(nSpikes,kernel,'same')/4;

% Number of bins for x and y
nBinsX = nBins(1);
if length(nBins) == 1,
	nBinsY = nBinsX;
	nBins(2) = nBins;
else
	nBinsY = nBins(2);
end

% Min and max speed and acceleration
xLims = [0 limits(1)];
yLims = limits(2:3);

% Bin rate and acceleration
rate = Bin(data.rate,xLims,nBinsX);
acceleration = Bin(data.acceleration,yLims,nBinsY);

% Compute phase map
map.phase = angle(Smooth(Accumulate([rate acceleration],exp(j*data.phase),[nBinsX nBinsY]),smooth))';
map.count = Smooth(Accumulate([rate acceleration],1,[nBinsX nBinsY]),smooth)';

% Bins
map.x = linspace(xLims(1),xLims(2),nBinsX);
map.y = linspace(yLims(1),yLims(2),nBinsY);


%  figure;
%  hold on;
%  PlotXY(phases);
%  PlotTicks(spikes,'size',1,'k');
%  plot(start,nSpikes,'b*');
%  plot(start,data.rate,'r*');
