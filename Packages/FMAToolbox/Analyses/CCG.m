%CCG - Compute multiple cross- and auto-correlograms
%
%  USAGE
%
%    [ccg,t] = CCG(times,groups,<options>)
%
%    times          times of all events (sorted)
%    groups         group IDs for each event in time list
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%    =========================================================================
%
%  SEE
%
%    See also ShortTimeCCG.

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [ccg,t] = CCG(times,groups,varargin)

% Default values
duration = 2;
binSize = 0.01;

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdvector(times),
	error('Parameter ''times'' is not a real-valued vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(groups) && ~isdvector(groups),
	error('Parameter ''groups'' is not a real-valued scalar or vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(groups) && length(times) ~= length(groups),
	error('Parameters ''times'' and ''groups'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
groups = groups(:);
times = times(:);

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
  end
end

% Number of groups, number of bins, etc.
if length(groups) == 1,
	groups = ones(length(times),1);
	nGroups = 1;
else
	nGroups = max(unique(groups));
end
halfBins = round(duration/binSize/2);
nBins = 2*halfBins+1;
t = (-halfBins:halfBins)'*binSize;

if length(times) <= 1,
	ccg = zeros(nBins,nGroups,nGroups);
	return
end

% Compute CCGs
nEvents = length(times);
counts = double(CCGHeart(times,uint32(groups),binSize,uint32(halfBins)));

% Reshape the results
n = max(groups);
counts = reshape(counts,[nBins n n]);
if n < nGroups,
	counts(nBins,nGroups,nGroups) = 0;
end

% Count events in each group
for g = 1:nGroups,
	nEventsPerGroup(g) = sum(groups==g);
end

% Reshape
for g1 = 1:nGroups,
	for g2 = g1:nGroups,
		ccg(:,g1,g2) = flipud(counts(:,g1,g2));
	end
end
