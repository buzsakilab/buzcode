function tuned = TuneArtefactTimes(x,times,varargin)

%TuneArtefactTimes - Fine-tune artefact (e.g. stimulation) times.
%
%  Fine tuning can be done using either the LFP (faster), or the wideband
%  data (slower but more accurate). The latter is especially useful in cases
%  where the artefacts are no longer visible in the LFP signal due to low-
%  pass filtering.
%
%  USAGE
%
%    % Using LFP (pass the actual samples)
%    times = TuneArtefactTimes(lfp,artefacts,<options>)
%
%    % Using wideband data (pass the channel ID because the data would not fit
%    % in memory, so it will have to be read piecewise from disk)
%    times = TuneArtefactTimes(channel,artefacts,<options>)
%
%    lfp            local field potentials
%    channel        recording channel (wideband data)
%    artefacts      approximate artefact times
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'method'      artefacts can be detected using either the derivative
%                   ('slope', default) or the amplitude ('amplitude') of the
%                   signal
%     'durations'   interval (in s) around each approximate artefact time
%                   where the actual artefact should be sought
%                   (default = [-0.0075 0.0075])
%     'verbose'     display information about ongoing processing
%                   (default = 'off')
%    =========================================================================
%
%  SEE
%
%    See also RemoveArtefacts.

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
durations = [-0.0075 0.0075];
verbose = 'off';
nChunks = 100;
method = 'slope';

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help TuneArtefactTimes">TuneArtefactTimes</a>'' for details).');
end

% LFP or wideband?
if isiscalar(x),
	source = 'wideband';
	channel = x;
else
	source = 'lfp';
	lfp = x;
	if size(lfp,2) ~= 2,
		error('Parameter ''lfp'' is not a Nx2 matrix (type ''help <a href="matlab:help TuneArtefactTimes">TuneArtefactTimes</a>'' for details).');
	end
end

% Check parameter sizes
if ~isdvector(times),
	error('Parameter ''times'' is not a vector (type ''help <a href="matlab:help TuneArtefactTimes">TuneArtefactTimes</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help TuneArtefactTimes">TuneArtefactTimes</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','>'),
			error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help TuneArtefactTimes">TuneArtefactTimes</a>'' for details).');
			end
		case 'method',
			method = lower(varargin{i+1});
			if ~isstring_FMAT(method,'slope','amplitude'),
				error('Incorrect value for property ''method'' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
			end
		case 'verbose',
			verbose = lower(varargin{i+1});
			if ~isstring_FMAT(verbose,'on','off'),
				error('Incorrect value for property ''verbose'' (type ''help <a href="matlab:help Sync">Sync</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help TuneArtefactTimes">TuneArtefactTimes</a>'' for details).']);

	end
end

% Discard duplicates
cleaned = unique(times);
if length(times) ~= length(cleaned),
	warning('Data contains duplicates - these will be discarded.');
end
times = cleaned;

if strcmp(source,'lfp'),
	% Process LFP all at once
	% Time intervals where the artefacts should be sought
	intervals = [times+durations(1) times+durations(2)];
	tuned = process(lfp,times,intervals,verbose);
else
	% Split data in chunks
	start = times(1) + durations(1);
	stop = times(end) + durations(2);
	totalTime = stop-start;
	dt = totalTime/nChunks;
	tuned = [];
	for i = 1:nChunks,
		% Time intervals where the artefacts should be sought
		t = Restrict(times,[start start+dt]);
		intervals = [t+durations(1) t+durations(2)];
		data = GetWidebandData(channel,'intervals',[start start+dt]);
		tuned = [tuned;process(data,t,intervals,method,verbose)];
		start = start + dt;
	end
end

% Actual processing

function tuned = process(signal,times,intervals,method,verbose)

% Differentiate signal
if strcmp(method,'slope'),
	signal(:,2) = [0;diff(signal(:,2))];
end

% Select the samples within these intervals
[in,i,j] = InIntervals(signal(:,1),intervals,'verbose',verbose);
% Transform this into a matrix where each line contains the data for one artefact, and find minima
m = full(sparse(i(i~=0),j(j~=0),signal(in,2)));
infs = ~logical(full(sparse(i(i~=0),j(j~=0),1)));
m(infs) = inf;
[unused,where] = min(m,[],2);
% Build the same matrix for data timestamps
t = full(sparse(i(i~=0),j(j~=0),signal(in,1)));
% Times
s1 = size(m,1);
s2 = size(m,2);
where = logical(full(sparse(1:s1,where,1,s1,s2)));
tuned = sort(t(where));
