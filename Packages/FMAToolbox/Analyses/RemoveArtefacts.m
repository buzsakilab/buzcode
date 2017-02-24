function cleaned = RemoveArtefacts(lfp,times,varargin)

%RemoveArtefacts - Remove artefacts from local field potentials.
%
% Remove large deflections such as stimulation artefacts from local field
% potentials, i.e. flatten the data around the times of artefacts.
%
%  USAGE
%
%    cleaned = RemoveArtefacts(lfp,times,<options>)
%
%    lfp            local field potential given as (t,V) pairs
%    times          timestamps of the artefacts
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'durations'   pre and post durations (in s) of the artefacts
%                   (default = [-0.001 0.002])
%    =========================================================================
%
%  SEE
%
%    See also TuneArtefactTimes.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
durations = 0.001 * [-1 2];

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help RemoveArtefacts">RemoveArtefacts</a>'' for details).');
end

% Check parameter sizes
if size(lfp,2) < 2,
	error('Parameter ''lfp'' is not a matrix (type ''help <a href="matlab:help RemoveArtefacts">RemoveArtefacts</a>'' for details).');
end
if size(times,2) ~= 1,
	error('Parameter ''times'' is not a vector (type ''help <a href="matlab:help RemoveArtefacts">RemoveArtefacts</a>'' for details).');
end

% Parse parameter list
for j = 1:2:length(varargin),
	if ~ischar(varargin{j}),
		error(['Parameter ' num2str(j+7) ' is not a property (type ''help <a href="matlab:help RemoveArtefacts">RemoveArtefacts</a>'' for details).']);
	end
	switch(lower(varargin{j})),
		case 'durations',
			durations = varargin{j+1};
			if ~isdvector(durations,'#2','<'),
				error('Incorrect value for ''durations'' (type ''help <a href="matlab:help RemoveArtefacts">RemoveArtefacts</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help RemoveArtefacts">RemoveArtefacts</a>'' for details).']);

	end
end

intervals = [times+durations(1) times+durations(2)];
in = InIntervals(lfp,intervals);
cleaned = [lfp(:,1) interp1(lfp(~in,1),lfp(~in,2),lfp(:,1))];
