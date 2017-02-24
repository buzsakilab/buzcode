function spikes = GetSpikeTimes(units,varargin)

%GetSpikeTimes - Get spike timestamps.
%
%  USAGE
%
%    spikes = GetSpikeTimes(units,<options>)
%
%    units          optional list of units, i.e. [electrode group, cluster] pairs;
%                   special conventions:
%                     cluster = -1   all clusters
%                     cluster = -2   all clusters except artefacts (cluster 0)
%                     cluster = -3   all clusters except artefacts (cluster 0)
%                                    and MUA (cluster 1)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'output'      'time' returns only timestamps, 'full' lists electrode
%                   group and cluster for each spike (default = 'time')
%    =========================================================================
%
%  EXAMPLES
%
%    % timestamps for all spikes
%    s = GetSpikeTimes;
%
%    % timestamps for units [1 7] and [4 3]
%    s = GetSpikeTimes([1 7;4 3]);
%
%    % timestamps for all units on electrode group 5 and unit [6 3]
%    s = GetSpikeTimes([5 -1;6 3]);
%
%    % timestamps for all units on electrode group 5, except artefacts
%    s = GetSpikeTimes([5 -2]);
%
%    % timestamps, electrode groups and clusters, for all spikes
%    s = GetSpikeTimes('output','full');
%
%  NOTE
%
%    An electrode group is an ensemble of closely spaced electrodes that record from
%    the same neurons, e.g. a single wire electrode, or a wire tetrode, or a multisite
%    silicon probe, etc.
%
%  SEE
%
%    See also LoadSpikeTimes, PlotTicks.


% Copyright (C) 2004-2011 by Michaël Zugaro
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
output = 'time';

% Optional parameter
if ischar(units),
	varargin = {units,varargin{:}};
	units = []; % all
else
	if ~isempty(units) && (~isimatrix(units) || size(units,2) ~= 2),
		error('Incorrect list of units (type ''help <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>'' for details).');
	end
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'output',
			output = lower(varargin{i+1});
			if ~isstring_FMAT(output,'time','full'),
				error('Incorrect value for property ''output'' (type ''help <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>'' for details).']);

  end
end

spikes = DATA.spikes;
if isempty(spikes), return; end

% Selected units only
if ~isempty(units),
	nUnits = size(units,1);
	selected = zeros(size(spikes(:,1)));
	for i = 1:nUnits,
		channel = units(i,1);
		cluster = units(i,2);
		switch cluster,
			case -1,
				selected = selected | spikes(:,2) == channel;
			case -2,
				selected = selected | (spikes(:,2) == channel & spikes(:,3) ~= 0);
			case -3,
				selected = selected | (spikes(:,2) == channel & spikes(:,3) ~= 0 & spikes(:,3) ~= 1);
			otherwise,
				selected = selected | (spikes(:,2) == channel & spikes(:,3) == cluster);
		end
	end
	spikes = spikes(selected,:);
end

if strcmp(output,'time'),
	spikes = spikes(:,1);
end
