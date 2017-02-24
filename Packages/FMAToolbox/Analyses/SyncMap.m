function [map,timeBins] = SyncMap(synchronized,indices,varargin)

%SyncMap - Create a map from successive event-synchronized data.
%
% Using N repetitions of similar data centered around synchronizing
% events (e.g. N evoked potentials), create a map v = f(t,i) where
% t is the time relative to the synchronizing events and i is the
% occurrence (from 1 to N).
%
%  USAGE
%
%    [map,timeBins] = SyncMap(synchronized,indices,<options>)
%
%    synchronized   event-synchronized samples
%    indices        list of synchronizing event indices for each sample
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'durations'   durations before and after synchronizing events for each
%                   trial (in s) (default = [-0.5 0.5])
%     'nBins'       total number of time bins (default 100)
%     'smooth'      smoothing size (0 = no smoothing) (default = 0.01*nBins)
%    =========================================================================
%
%  SEE
%
%    See also Sync, SyncHist, PlotSync, PETHTransition.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
durations = [-0.5 0.5];
nBins = 100;
smooth = [];

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).');
end

% Check parameter sizes
if size(synchronized,2) > 2,
	error('Parameter ''synchronized'' is not a Nx2 matrix (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).');
end
if size(indices,2) ~= 1,
	error('Parameter ''indices'' is not a vector (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).');
end
if size(indices,1) ~= size(synchronized,1),
	error('Parameters ''synchronized'' and ''indices'' have different lengths (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','<'),
				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).');
			end
		case 'nbins',
			nBins = varargin{i+1};
			if ~isiscalar(nBins,'>0'),
				error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).');
			end
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdvector(smooth,'>=0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SyncMap">SyncMap</a>'' for details).']);
	end
end


if isempty(smooth),
	smooth = ceil(0.01*nBins);
end

start = durations(1);
stop = durations(2);
timeBinSize = (stop - start)/nBins;
timeBins = (start:timeBinSize:stop-timeBinSize)+timeBinSize/2;
binnedTime = Bin(synchronized(:,1),durations,nBins);
if size(synchronized,2) == 1,
	% Occurrences
	s = Accumulate([indices binnedTime],1);
	map = Smooth(s,smooth);
else
	% Values
	s = Accumulate([indices binnedTime],synchronized(:,2));
	n = Accumulate([indices binnedTime],1);
	n(n==0) = 1;
	map = Smooth(s,[smooth 0])./Smooth(n,[smooth 0]);
end
