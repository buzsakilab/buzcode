function stats = ReconstructPosition(positions,spikes,phases,varargin)

%ReconstructPosition - Bayesian reconstruction of positions from spike trains.
%
% Instantaneous positions are reconstructed using a Bayesian algorithm.
% Instantaneous population firing rates can be estimated either over fixed time
% windows, or over fractions of the theta cycle (or of any other brain rhythm).
% Similarly, positions will be reconstructed either over time or phase windows.
% The model is first trained on a subset of the data, then tested on the rest.
%
% USAGE
%
%    stats = ReconstructPosition(positions,spikes,phases,<options>)
%
%    positions      linear or two-dimensional positions, in [0..1]
%    spikes         list of (t,group,cluster) triplets (obtained via e.g.
%                   <a href="matlab:help GetSpikes">GetSpikes</a>, using full output)
%    phases         optional unwrapped phase of the LFP (see <a href="matlab:help Phase">Phase</a>)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'training'    time interval over which the model should be trained
%                   (default = first half of the position data)
%     'window'      length of the time or phase window (default = 0.020 s for
%                   time, and pi/3 for phases)
%     'type'        two letters (one for X and one for Y) indicating which
%                   coordinates are linear ('l') and which are circular ('c')
%                   - for 1D data, only one letter is used (default 'll')
%     'nBins'       firing curve or map resolution (default = [200 200])
%    =========================================================================
%
%   OUTPUT
%
%     stats.positions     real position across time or phase windows
%     stats.spikes        cell firing vector across time or phase windows
%     stats.estimations   estimated position across time or phase windows
%     stats.errors        estimation error across time or phase windows
%     stats.average       average estimation error in each phase window
%     stats.windows       time windows (possibly computed from phases)
%     stats.phases        phase windows (empty for fixed time windows)
%

% Copyright (C) 2012-2013 by Michaël Zugaro, (C) 2012 by Karim El Kanbi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
wt = 0.020; % default time window
wp = pi/3; % default phase window
window = [];
nBins = 200;
training = 0.5;
type = '';
nDimensions = 1;

% Check number of parameters
if nargin < 2 || mod(length(varargin),2) ~= 0,
	builtin('error','Incorrect number of parameters (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end

% Optional parameter 'phases'
if nargin == 2,
	phases = [];
elseif nargin >= 3 && ischar(phases),
	varargin = {phases,varargin{:}};
	phases = [];
end

% Check parameter sizes
if ~isdmatrix(positions),
	builtin('error','Incorrect positions (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end
if ~isdmatrix(spikes),
	builtin('error','Incorrect spikes (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end
if size(positions,2) > 3,
	nDimensions = 2;
end
if ~isempty(phases) && ~isdmatrix(phases),
	builtin('error','Incorrect value for property ''phases'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end

% Check number of output parameters
if isempty(phases) && nargout == 3,
	builtin('error','Too many output parameters or missing phases (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end

% Parse parameters
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		builtin('error',['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'training',
			training = varargin{i+1};
			if ~isdvector(training,'>'),
				builtin('error','Incorrect value for property ''training'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				builtin('error','Incorrect value for property ''window'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				builtin('error','Incorrect value for property ''show'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'nBins',
			nBins = varargin{i+1};
			if ~isinteger(nBins),
				builtin('error','Incorrect value for property ''nBins'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'type',
			type = lower(varargin{i+1});
			if (nDimensions == 1 && ~isstring_FMAT(type,'cc','cl','lc','ll')) || (nDimensions == 2 && ~isstring_FMAT(type,'ccl','cll','lcl','lll','ccc','clc','lcc','llc')),
				builtin('error','Incorrect value for property ''type'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		otherwise,
			builtin('error',['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
	end
end

% Defaults
if isempty(window),
	if isempty(phases),
		window = wt;
	else
		window = wp;
		if ~isiscalar((2*pi)/window),
			builtin('error',['Incorrect phase window: not an integer fraction of 2pi (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
		end
	end
end
if isdscalar(training),
	training = [-Inf positions(1,1)+training*(positions(end,1)-positions(1,1))];
else
	if min([positions(1,1) spikes(1,1)]) < training(1) & min([positions(end,1) spikes(end,1)]) > training(2),
		builtin('error',['Spikes or positions occur both before and after training interval (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
	end
end
if isempty(type),
	if nDimensions == 2,
		type = 'lll';
	else
		type = 'll';
	end
end
nBinsX = nBins(1);
if length(nBins) > 2,
	nBinsY = nBins(2);
else
	if nDimensions == 2,
		nBinsY = nBinsX;
	else
		nBinsY = 1;
	end
end


% List units, assign them an ID (number them from 1 to N), and associate these IDs with each spike
% (IDs will be easier to manipulate than (group,cluster) pairs in subsequent computations)
[units,~,i] = unique(spikes(:,2:end),'rows');
nUnits = length(units);
index = 1:nUnits;
id = index(i)';
spikes = [spikes(:,1) id];

% Split data (training vs test)
trainingPositions = Restrict(positions,training);
trainingSpikes = Restrict(spikes,training);
testPositions = positions(~InIntervals(positions,training),:);
testSpikes = spikes(~InIntervals(spikes,training),:);

% TRAINING

% Compute occupancy probability P(x) (i.e. normalized occupancy map)
firstUnit = trainingSpikes(:,2) == 1;
s = trainingSpikes(firstUnit,1);
map = Map(trainingPositions,s,'nbins',nBins,'smooth',5,'type',type);
Px = map.time;
Px = Px ./ sum(Px(:));

% Compute average firing probability lambda for each unit (i.e. firing maps)
lambda(:,:,1) = map.z;
for i = 2:nUnits,
	nextUnit = trainingSpikes(:,2) == i;
	s = trainingSpikes(nextUnit,1);
	map = Map(trainingPositions,s,'nbins',nBins,'smooth',5,'type',type);
	lambda(:,:,i) = map.z;
end

% TEST

% Determine time windows (using unwrapped phases if necessary)
if ~isempty(phases),
	testPhases = phases(~InIntervals(phases,training),:);
	drop = testPhases(:,1) < testPositions(1,1);
	testPhases(drop,:) = [];
	startPhase = ceil(testPhases(1,2)/(2*pi))*2*pi;
	stopPhase = floor(testPhases(end,2)/(2*pi))*2*pi;
	windows = (startPhase:window:stopPhase)';
	stats.phases = windows;
	windows = Interpolate(testPhases(:,[2 1]),windows);
	windows = [windows(1:end-1,2) windows(2:end,2)];
else
	stats.phases = [];
	windows = (testPositions(1,1):window:testPositions(end,1))';
	windows = [windows(1:end-1) windows(2:end)];
end
nWindows = size(windows,1);

stats.estimations = nan(nBinsY,nBinsX,nWindows);
stats.spikes = zeros(nUnits,nWindows);
% Loop over data windows
for i = 1:nWindows,
    
	% Get spikes for this window
	s = Restrict(testSpikes,windows(i,:));

	if isempty(s),
		% No spikes: set uniform probability
		stats.estimations(:,:,i) = ones(nBinsY,nBinsX,1)/(nBinsX*nBinsY);
		continue;
	end
    
	% Population spike count vector
	stats.spikes(:,i) = Accumulate(s(:,2),1,nUnits);
	% To avoid 'for' loops, prepare for vector computation:
	% assign a spike count to each position and unit (3D array)
	n = reshape(repmat(stats.spikes(:,i),1,nBinsX*nBinsY)',nBinsY,nBinsX,nUnits);

	% For each cell i, compute P(ni|x) using a Poisson model
	dt = windows(i,2) - windows(i,1);
	Pnix = (dt*lambda).^n./factorial(n).*exp(-dt*lambda);
	% Compute P(n|x) assuming independent probabilities across units (hmm...)
	% i.e. P(n|x) = product over i of P(ni|x)
	Pnx = prod(Pnix,3);
            
	% Compute P(n) = sum over x of P(n|x)*P(x)
	Pn = sum(sum(Pnx.*Px));
        
	% Compute P(x|n) = P(n|x)*P(x)/P(n)
	Pxn = Pnx .* Px / Pn;

	% Store result
	stats.estimations(:,:,i) = Pxn;

end
stats.estimations = squeeze(stats.estimations);
stats.windows = windows;

% Estimation error

stats.errors = [];
stats.average = [];
if nDimensions == 1,
	% Bin test positions and compute distance to center
	stats.positions = Interpolate(testPositions,windows(:,1));
	stats.positions(:,2) = Bin(stats.positions(:,2),[0 1],nBinsX);
	dx = (round(nBinsX/2)-stats.positions(:,2))';
	% Shift estimated position by the real distance to center
	stats.errors = CircularShift(stats.estimations(:,1:length(dx)),dx);
	% Average over one or more cycles
	k = 2*pi/window;
	n = floor(size(stats.errors,2)/k)*k;
	stats.average = reshape(stats.errors(:,1:n),nBins,k,[]);
	stats.average = nanmean(stats.average,3);
else
	warning('Computation of estimation error not yet implemented for 2D environments');
end
