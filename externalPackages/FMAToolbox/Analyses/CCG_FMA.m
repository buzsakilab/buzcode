function [ccg,t,tau,c] = CCG_FMA(times,id,varargin)

%CCG - Compute multiple cross- and auto-correlograms or cross-covariances
%
%  USAGE
%
%    [ccg,t,tau,c] = CCG(times,id,<options>)
%
%    times          times of all events (see NOTE below)
%    id             ID for each event (e.g. unit ID) from 1 to n
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%     'smooth'      smoothing size in bins (0 = no smoothing, default)
%     'groups'      group number (1 or 2) for each event, used to restrict
%                   cross-correlograms to pairs across two groups of events
%                   (see EXAMPLE #2 below)
%     'mode'        'ccg' or 'ccv' (default = 'ccg')
%     'alpha'       significance level to determine correlated pairs
%     'totaltime'   total recording duration in s (default suppose only one 
%                   continuous recording block = max(times) - min(times))
%    =========================================================================
%
%  OUTPUT
%      ccg          value of cross-correlograms or cross-covariances
%                   dimensions are (nbins,m,n) where m is the number of
%                   reference time series (e.g. reference units) and n the
%                   number of referenced time series (in general m = n,
%                   except when using option 'groups')
%      t            time bins
%      tau          lag times for a each pair (mode 'ccv' only)
%      c            maximum cross-covariance for a each pair (mode 'ccv' only)
%
%
%  NOTE
%
%    Parameters 'times', 'id' and 'group' can be obtained using <a href="matlab:help CCGParameters">CCGParameters</a>.
%    As a special case, when computing the correlograms of spike trains, one
%    can use the output of <a href="matlab:help GetSpikes">GetSpikes</a> either directly or in combination with
%    <a href="matlab:help CCGParameters">CCGParameters</a>. See EXAMPLES below.
%
%  EXAMPLES
%
%    % Auto- and cross-correlograms between all neurons
%    spikes = GetSpikes('output','numbered');
%    [ccg,t] = CCG(spikes(:,1),spikes(:,2));
%
%    % Only tetrode #1 vs tetrode #2 (e.g. mPFC vs HPC neurons)
%    pfc = GetSpikes([1 -1],'output','numbered');
%    hpc = GetSpikes([2 -1],'output','numbered');
%    [s,ids,groups] = CCGParameters(pfc,hpc,2);
%    [ccg,t] = CCG(s,ids,'groups',groups);
%
%    % Between stimulations and MUA spikes
%    spikes = GetSpikes;
%    stimulatios = GetEvents('Stimulation');
%    d = [spikes(:,1) ones(size(spikes,1)) ; stimulations 2*ones(size(stimulations,1))];
%    d = sortrows(d);
%    [ccg,t] = CCG(d(:,1),d(:,2));
%
%    % To compute cross-covariances
%    [ccv,t,tau,C] = CCG(times,ids,'mode','ccv');
%
%  SEE
%
%    See also CCGParameters, ShortTimeCCG.

% Copyright (C) 2012-2013 by Michaël Zugaro, Marie Goutierre
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Default values
duration = 2;
binSize = 0.01;
smooth = 0;
groups = [];
mode = 'ccg';
alpha = 0.05;
totalTime = max(times)-min(times);

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdvector(times),
	error('Parameter ''times'' is not a real-valued vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(id) && ~isdvector(id),
	error('Parameter ''id'' is not a real-valued scalar or vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(id) && length(times) ~= length(id),
	error('Parameters ''times'' and ''id'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
id = id(:);
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
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>=0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'groups',
			groups = varargin{i+1};
			if ~isempty(groups) && ~isdvector(groups) && length(times) ~= length(groups)
				error('Incorrect value for property ''groups'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'alpha',
			alpha = varargin{i+1};
			if ~isdscalar(alpha,'>0'),
				error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'mode',
			mode = varargin{i+1};
			if ~isstring(mode,'ccg','ccv'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'totaltime',
			totalTime = varargin{i+1};
			if ~isdscalar(totalTime,'>0'),
				error('Incorrect value for property ''totaltime'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
end

tau = [];
c = [];

% Number of IDs, number of bins, etc.
if length(id) == 1,
	id = ones(length(times),1);
	nIDs = 1;
else
	nIDs = length(unique(id));
end
if nIDs ~= max(id),
	error('Incorrect IDs (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
halfBins = round(duration/binSize/2);
nBins = 2*halfBins+1;
t = (-halfBins:halfBins)'*binSize;

if length(times) <= 1,
	return
end

% Sort events in time and compute CCGs
[times,i] = sort(times);
id = id(i);
if ~isempty(groups),
	groups = groups(i);
end
counts = CCGEngine(times,id,binSize,halfBins);

% Reshape the results
n = max(id);
counts = reshape(counts,[nBins n n]);
if n < nIDs,
	counts(nBins,nIDs,nIDs) = 0;
end

% Restrict the results to inter-group CCGs if requested
if ~isempty(groups),
	group1 = unique(id(groups == 1));
	group2 = unique(id(groups == 2));
	nGroup1 = length(group1);
	nGroup2 = length(group2);
	ccg = zeros(nBins,nGroup1,nGroup2);
	for i = 1:nGroup1,
		for j = 1:nGroup2,
			ccg(:,i,j) = Smooth(flipud(counts(:,group1(i),group2(j))),smooth);
		end
	end
else
	ccg = zeros(nBins,nIDs,nIDs);
	% Compute corr(A,B) for each unique unordered pair (A,B)
	for g1 = 1:nIDs,
		for g2 = g1:nIDs,
			ccg(:,g1,g2) = Smooth(flipud(counts(:,g1,g2)),smooth);
		end
	end
	% corr(B,A) and corr(B,A) symmetric
	for g1 = 1:nIDs,
		for g2 = 1:g1-1,
			ccg(:,g1,g2) = flipud(squeeze(ccg(:,g2,g1)));
		end
	end
end


if strcmp(mode,'ccv'),

	% Determine mean event rate for each ID
	eventRate = zeros(nIDs,1);
	for i = 1:nIDs,
		eventRate(i) = sum(id==i)/totalTime;
	end

	% Determine standardized cross-covariances
	ccv = zeros(size(ccg));
	tau = zeros(size(ccg,2),size(ccg,3));
	c = zeros(size(ccg,2),size(ccg,3));

	nPairs = size(ccg,2)*size(ccg,3);
	disp(['# pairs: ' int2str(nPairs)]);

	threshold = sqrt(2)*erfinv(1-(alpha/length(t)));

	for i = 1:size(ccg,2),
		for j = 1:size(ccg,3),
		
			% Compute and normalize CCVs from CCGs
			if ~isempty(groups),
				rate = eventRate(group1(i))*eventRate(group2(j));
			else
				rate = eventRate(i)*eventRate(j);
			end
			ccv(:,i,j) = sqrt((binSize*totalTime)/rate) * (ccg(:,i,j)/(binSize*totalTime)-rate);

			% Smooth with a 3-bin boxcar
			data = ccv(:,i,j);
			top = flipud(data(1:size(ccg,1),:));
			bottom = flipud(data(end-size(ccg,1)+1:end,:));
			data = [top;data;bottom];
			data = filter([1 1 1],3,data);
			n = size(data,1);
			d = n - size(ccg,1);
			start = d/2+1;
			stop = start + size(ccg,1) - 1;
			ccv(:,i,j) = data(start:stop);

			% Find the peak lag time and value
			[~,maxIndex] = max(ccv(:,i,j));
			tau(i,j) = median(t(maxIndex));
			c(i,j) = median(ccv(maxIndex,i,j));
			% Previous version of the code (discard?)
			% [~,maxIndex] = max(ccv(:,i,j));
			% tau(i,j) = t(maxIndex);
			% c(i,j) = median(ccv(max(1,maxIndex-3):min(end,maxIndex+3),i,j));

			% Keep only significantly correlated pairs
			if ~any(abs(ccv(:,i,j))>threshold),
				tau(i,j) = NaN;
			end
			
		end
	end

	nCorrelatedPairs = sum(~isnan(tau(:)));
	disp(['# significantly correlated pairs: ' int2str(nCorrelatedPairs)]);

	ccg = ccv;
	
end




%    % Only tetrode #1 vs tetrode #2 (e.g. mPFC vs HPC neurons)
%    pfc = GetSpikes([1 -1],'output','numbered');
%    hpc = GetSpikes([2 -1],'output','numbered');
%    m = max(pfc(:,2));
%    [spikes,i] = sortrows([pfc(:,1);hpc(:,1)]);
%    ids = [pfc(:,2);hpc(:,2)+m];
%    ids = ids(i);
%    groups = [ones(size(pfc(:,1)));2*ones(size(hpc(:,1)))];
%    groups = groups(i);
%    [ccg,t] = CCG(s,ids,'groups',groups);
