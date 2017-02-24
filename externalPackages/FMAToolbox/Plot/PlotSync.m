function PlotSync(synchronized,indices,varargin)

%PlotSync - Plot successive occurrences of a multidimensional variable.
%
% Plot resynchronized samples such as spike rasters or successive evoked
% potentials.
%
%  USAGE
%
%    PlotSync(synchronized,indices,<options>)
%
%    synchronized   synchronized samples (obtained using <a href="matlab:help Sync">Sync</a>)
%    indices        synchronizing event indices (obtained using <a href="matlab:help Sync">Sync</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'durations'   durations before and after synchronizing events for each
%                   trial (in s) (default = [-0.5 0.5])
%     'spacing'     vertical spacing between samples occurring around
%                   successive synchronizing events (default = 1 for point
%                   processes, 0 for continuous data)
%    =========================================================================
%
%  SEE
%
%    See also Sync, SyncHist, SyncMap.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check arguments
if nargin < 2 | mod(length(varargin),2) ~= 0,
error('Incorrect number of parameters (type ''help <a href="matlab:help PlotSync">PlotSync</a>'' for details).');
end
if isempty(indices) || isempty(synchronized), return; end
if min(size(indices)) ~= 1,
	error('''indices'' is not a vector (type ''help <a href="matlab:help PlotSync">PlotSync</a>'' for details).');
end
if length(synchronized) ~= length(indices),
	error('''synchronized'' and ''indices'' have different lengths (type ''help <a href="matlab:help PlotSync">PlotSync</a>'' for details).');
end

% Default values
durations = [-0.5 0.5];
if size(synchronized,2) == 1,
	% Point process data
	pointProcess = true;
	spacing = 1;
	synchronized(:,2) = 1;
else
	% Continuous-valued data
	pointProcess = false;
	spacing = 0;
end


% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotSync">PlotSync</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','<'),
				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help PlotSync">PlotSync</a>'' for details).');
			end

		case 'spacing',
			spacing = varargin{i+1};
			if ~isdscalar(spacing),
				error('Incorrect value for property ''spacing'' (type ''help <a href="matlab:help PlotSync">PlotSync</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotSync">PlotSync</a>'' for details).']);
	end
end

% Plot
hold on;
if pointProcess,
	PlotTicks(synchronized(:,1),(indices-1)*spacing);
else
	nSync = max(indices);
	for i = 1:nSync,
		trial = indices==i;
		p = plot(synchronized(trial,1),(i-1)*spacing+synchronized(trial,2));
	end
end
% Determine axes limits
minY = min(synchronized(:,2)+spacing*(indices-1));
maxY = max(synchronized(:,2)+spacing*(indices-1));
if minY == maxY, minY = minY-eps; end
if min(synchronized(:,1)) < durations(1) | max(synchronized(:,1)) > durations(2),
	warning(['Some data points do not fall within [' num2str(durations(1)) ';' num2str(durations(2)) ']']);
end

% Update axes
set(gca,'XLim',durations,'YLim',[minY maxY]);
line([0 0],ylim,'linestyle','--','color',[0 0 0]);
