function spikes = ThresholdSpikes(amplitudes,factor,varargin)

%ThresholdSpikes - Post-hoc threshold correction for spike detection.
%
%  USAGE
%
%    spikes = ThresholdSpikes(amplitudes,factor,<options>)
%
%    amplitudes     [time, electrode group, cluster, amplitude] tuples
%                   e.g. obtained using <a href="matlab:help GetSpikeAmplitudes">GetSpikeAmplitudes</a>
%    factor         the new threshold will be computed as the former threshold
%                   multiplied by this factor
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'units'       existing units (this is used to generate new cluster IDs)
%    =========================================================================
%
%  OUTPUT
%
%    The subset of spikes exceeding the new threshold. The cluster IDs are
%    unchanged, unless the optional parameter 'units' is provided, in which
%    case each unit is assigned a new cluster ID.
%
%  SEE
%
%    See also GetSpikeAmplitudes.

% Copyright (C) 2009-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
allUnits = [];

% Parse parameters
if nargin < 2,
	error('Missing parameters (type ''help <a href="matlab:help ThresholdSpikes">ThresholdSpikes</a>'' for details).');
end
if ~isdmatrix(amplitudes),
	error('Incorrect amplitudes (type ''help <a href="matlab:help ThresholdSpikes">ThresholdSpikes</a>'' for details).');
end
if ~isdscalar(factor),
	error('Incorrect factor (type ''help <a href="matlab:help ThresholdSpikes">ThresholdSpikes</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ThresholdSpikes">ThresholdSpikes</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help ThresholdSpikes">ThresholdSpikes</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'units',
			allUnits = varargin{i+1};
			if ~isdmatrix(allUnits,'>=0') || size(allUnits,2) ~= 2,
				error('Incorrect value for property ''units'' (type ''help <a href="matlab:help ThresholdSpikes">ThresholdSpikes</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ThresholdSpikes">ThresholdSpikes</a>'' for details).']);
	end
end

% Which units are we processing?
units = unique(amplitudes(:,2:3),'rows');

spikes = [];
for i = 1:size(units,1),
	group = units(i,1);
	cluster = units(i,2);
	% Amplitudes for this unit
	this = amplitudes(amplitudes(:,2)==group&amplitudes(:,3)==cluster,:);
	% Guess former threshold
	peak = max(abs(this(:,4:end)),[],2);
	threshold = min(peak);
	% Update
	threshold = factor*threshold;
	% Spikes to keep
	keep = this(peak>=threshold,1:3);
	if ~isempty(allUnits),
		% New cluster ID = largest existing ID + 1
		newCluster = max(allUnits(allUnits(:,1)==group,2)) + 1;
		allUnits = [allUnits;group newCluster];
		keep(:,3) = newCluster;
		disp(['Unit (' int2str(group) ',' int2str(cluster) ') -> (' int2str(group) ',' int2str(newCluster) ')']);
	end
	spikes = [spikes;keep];
end
