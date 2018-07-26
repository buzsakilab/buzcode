function [synchronized,indices] = Sync(samples,sync,varargin)

%Sync - Make sample timestamps relative to synchronizing events.
%
% Select samples that fall around synchronizing events, and make their
% timestamps relative to the synchronizing events. This can be used to
% build e.g. spike raster plots or successive evoked potentials.
%
%  USAGE
%
%    [synchronized,indices] = Sync(samples,sync,<options>)
%
%    samples        either a vector of timestamps (for a point process) or
%                   an Mx(N+1) matrix containing M (timestamp,values)
%                   (N+1)-tuples (for N time-varying measures).
%    sync           timestamps to synchronize on (e.g., brain stimulations)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'durations'   durations before and after synchronizing events for each
%                   trial (in s) (default = [-0.5 0.5])
%     'verbose'     display information about ongoing processing
%                   (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    synchronized   resynchronized samples
%    indices        for each sample, the corresponding synchronizing event index
%
%  EXAMPLE
%
%    [raster,indices] = Sync(spikes,stimuli);     % compute spike raster data
%    figure;PlotSync(raster,indices);             % plot spike raster
%
%  SEE
%
%    See also SyncHist, SyncMap, PlotSync, PETHTransition.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
durations = [-.5 .5];
verbose = false;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
end

% Check parameter sizes
if ~isdvector(sync),
	error('Parameter ''sync'' is not a vector (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
end
if isempty(sync) || isempty(samples),
	synchronized = [];
	indices = [];
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Sync">Sync</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','<='),
				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
			end
		case 'verbose',
			verbose = varargin{i+1};
			if ~isstring_FMAT(verbose,'on','off'),
				error('Incorrect value for property ''verbose'' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
			end
			verbose = strcmp(verbose,'on');
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).']);
	end
end

nSync = length(sync);

% disp([num2str(nSync) ' synchronizing events to process...']) ;

% Output matrices will be allocated in arbitrarily large blocks (and will be trimmed down later)
% This will be much faster than increasing matrix size by the exact appropriate amount for each synchronizing event
blockLength = 1e6;
synchronized = nan(blockLength,size(samples,2));
indices = nan(blockLength,1);

% This variable is used to speed up computations: to find points within each sync time window, start from the beginning
% of the previous sync window rather than from the very beginning of the sample data
previous = 1;

k = 1;
for i = 1:nSync,
	% Find samples within time window around this synchronizing event
	j = FindInInterval(samples,[sync(i)+durations(1) sync(i)+durations(2)],previous);
	if ~isempty(j),
		previous = j(1);
		j = (j(1):j(2))';
		nj = length(j);
		if verbose, disp([' sync ' int2str(i) ' (t=' num2str(sync(i)) '): ' int2str(length(j)) ' samples']); end
		% Synchronize samples
		if ~isempty(j),
			if length(indices) < k + nj,
				% Enlarge matrices if necessary
				synchronized = [synchronized;nan(blockLength,size(samples,2))];
				indices = [indices;nan(blockLength,1)];
			end
			synchronized(k:k+nj-1,:) = [samples(j,1)-sync(i) samples(j,2:end)];
			indices(k:k+nj-1,1) = i*ones(size(j));
		end
		k = k + nj;
	end
end
indices(k:end) = [];
synchronized(k:end,:) = [];

