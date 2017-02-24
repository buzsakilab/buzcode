function amplitudes = GetSpikeAmplitudes(units,varargin)

%GetSpikeAmplitudes - Get spike amplitudes.
%
%  USAGE
%
%    amplitudes = GetSpikeAmplitudes(units,<options>)
%
%    units          list of [electrode group,cluster] pairs
%                   special conventions:
%                     cluster = -1   all clusters except artefacts and MUA
%                     cluster = -2   all clusters except artefacts
%                     cluster = -3   all clusters
%                   (artefacts are assumed to be in cluster 0, and MUA in 1)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'restrict'    list of time intervals to read from file
%    =========================================================================
%
%  SEE
%
%    See also GetSpikeTimes, GetSpikeWaveforms.

% Copyright (C) 2004-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Default values
intervals = [];

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help GetSpikeAmplitudes">GetSpikeAmplitudes</a>'' for details).');
end

if ~isimatrix(units,'>=-3') || size(units,2) ~= 2,
	error('Incorrect unit list (type ''help <a href="matlab:help GetSpikeAmplitudes">GetSpikeAmplitudes</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help GetSpikeAmplitudes">GetSpikeAmplitudes</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'restrict',
		intervals = varargin{i+1};
		if ~isdmatrix(intervals) || size(intervals,2) ~= 2,
			error('Incorrect value for property ''restrict'' (type ''help <a href="matlab:help GetSpikeAmplitudes">GetSpikeAmplitudes</a>'' for details).');
		end
		otherwise,
		error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GetSpikeAmplitudes">GetSpikeAmplitudes</a>'' for details).']);
	end
end

groups = unique(units(:,1));

% Check for incompatible list items
for i = 1:length(groups),
	group = groups(i);
	clusters = units(units(:,1)==group,2);
	if	(ismember(-1,clusters) && ~all(clusters==-1)) || (ismember(-2,clusters) && ~all(clusters==-2)) || (ismember(-3,clusters) && ~all(clusters==-3)),
		error('Incompatible list items (type ''help <a href="matlab:help GetSpikeAmplitudes">GetSpikeAmplitudes</a>'' for details).');
	end
end

amplitudes = [];
for i = 1:length(groups),

	% Clusters for this group
	group = groups(i);
	clusters = units(units(:,1)==group,2);
	
	% Load all waveforms for this group
	filename = [DATA.session.path '/' DATA.session.basename '.spk.' int2str(group)];
	nChannels = length(DATA.spikeGroups.groups{group});
	nSamplesPerWaveform = DATA.spikeGroups.nSamples(group);
	peak = DATA.spikeGroups.peakSamples(group);
	a = LoadSpikeAmplitudes(filename,nChannels,nSamplesPerWaveform,peak,DATA.rates.wideband);

	% Select appropriate clusters
	if ismember(-1,clusters),
		keep = a(:,3) ~= 0 & a(:,3) ~= 1;
	elseif ismember(-2,clusters),
		keep = a(:,3) ~= 0;
	elseif ismember(-3,clusters),
		keep = logical(ones(size(a,1),1));
	else
		keep = ismember(a(:,3),clusters);
	end
	a = a(keep,:);

	% Select timeframes
	if ~isempty(intervals),
		a = Restrict(a,intervals);
	end
	
	% Add to growing list
	amplitudes = [amplitudes;a];
	
end

% Sort list by time
amplitudes = sortrows(amplitudes);