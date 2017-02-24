function [ccg,x,y] = ShortTimeCCG(times1,times2,varargin)

%ShortTimeCCG - Time-varying auto/cross-correlograms of point processes.
%
% Determine how xcorrelograms vary in time, i.e. repeatedly compute xcorrelogram
% over successive short time windows.
%
%  USAGE
%
%    [ccg,bins,time] = ShortTimeCCG(times1,times2,<options>)
%
%    times1         first (reference) list of timestamps (in s)
%    times2         optional second list of timestamps (in s)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration (in s) of each xcorrelogram (default = 2)
%     'window'      duration (in s) of the time window over which each
%                   xcorrelogram is computed (default = 5*60)
%     'overlap'     overlap between successive windows (default = 0.8*window)
%     'smooth'      standard deviation for Gaussian kernel (default 0, no
%                   smoothing)
%     'mode'        'counts' yields raw event counts (default), and 'norm'
%                   normalizes each xcorrelogram to yield a probability
%                   distribution
%     'min'         discard time windows with fewer events than this threshold
%                   (default = 1)
%    =========================================================================
%
%  OUTPUT
%
%    ccg            MxN matrix, where each column is a xcorrelogram and each line
%                   is a time bin
%    bins           time bins for xcorrelograms
%    time           time
%
%  SEE
%
%    See also CCG, PlotShortTimeCCG.

% Copyright (C) 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% cross/auto correlogram?
auto = 0;

% Check parameter sizes
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
end
if nargin >= 2 && isa(times2,'char'),
	% Only one point process => autocorrelogram
	varargin = {times2,varargin{:}};
	times2 = times1;
	auto = 1;
end
if ~isvector(times1) || ~isvector(times2) || isempty(times1) || isempty(times2),
	error('Parameters ''times1'' and ''time2'' must be non-empty vectors (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
end

% Default values
binSize = 0.01;
duration = 2;
window = 5*60;
overlap = [];
mode = 'count';
smooth = [];
minEvents = 1;

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),

		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
			end

		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
			end

		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
			end

		case 'overlap',
			overlap = varargin{i+1};
			if ~isdscalar(overlap,'>0'),
				error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
			end

		case 'mode',
			mode = lower(varargin{i+1});
			if ~isstring_FMAT(mode,'norm','count'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
			end

		case 'smooth',
			smooth = varargin{i+1};
			if ~isdvector(smooth,'>=0') | length(smooth) > 2,
			error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
			end

		case 'min',
			minEvents = varargin{i+1};
			if ~isiscalar(minEvents,'>=0'),
				error('Incorrect value for property ''min'' (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>'' for details).']);

  end
end

% Default overlap?
if isempty(overlap), overlap = 0.8*window; end

% Get ready...
start = min([times1(1) times2(1)]+window/2);
stop = max([times1(end) times2(end)]-window/2);
times1 = times1(:);
times2 = times2(:);
times = [times1;times2];
groups = 1+[zeros(size(times1));ones(size(times2))];
halfBins = round(duration/binSize/2);

% Loop through data
i = 1;
t = start;
ccg = [];
while t+window/2 <= stop,
	x(i,1) = t;
	ok = InIntervals(times,t+[-0.5 0.5]*window);
	if sum(ok) < minEvents,
		ccg(1:(2*halfBins+1),i) = nan;
	else
		out = CCG(times(ok),groups(ok),'binSize',binSize,'duration',duration);
		if auto,
			ccg(:,i) = out(:,1,1);
		else
			ccg(:,i) = out(:,1,2);
		end
	end
	t = t + window - overlap;
	i = i + 1;
end

y = (-duration/2:binSize:duration/2)';

% Remove center bin for autocorrelograms
if auto,
	center = ceil(size(ccg,1)/2);
	ccg(center,:) = 0;
end

% Normalize?
if strcmp(mode,'norm'),
	ccg = ccg ./ repmat(sum(ccg,1),size(ccg,1),1);
end

% Smooth?
if ~isempty(smooth),
	ccg = Smooth(ccg,smooth);
end