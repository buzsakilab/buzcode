function [distance,units] = IsolationDistance(features,varargin)

%IsolationDistance - Determine the isolation quality for one or more clusters.
%
% Compute isolation distance, as defined by K. Harris (see Schmitzer-Torbert
% et al. 2005).
%
%  USAGE
%
%    [d,units] = IsolationDistance(features,<options>)
%
%    features       features obtained using e.g. <a href="matlab:help GetSpikeFeatures">GetSpikeFeatures</a> (or a
%                   subset thereof)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'units'       list of units to process (default = all clusters except 0
%                   and 1)
%     'show'        plot Mahalanobis distance PDF and CDF (default = 'off')
%    =========================================================================
%
%  OUTPUT
% 
%    d              For each cluster, a measure of the separation between the
%                   cluster and all other spikes. Namely, the Mahalanobis
%                   distance of the Nth closest spike outside the cluster
%                   (where N is the number of spikes inside the cluster).
%    units          corresponding (group,cluster) pairs
%

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
distance = [];
nFeatures = size(features,2)-3;
units = [];
show = 'off';

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help IsolationDistance">IsolationDistance</a>'' for details).');
end

% Check parameter sizes
if length(size(features)) ~= 2 | size(features,2) < 4,
	error('Incorrect parameter ''features'' (type ''help <a href="matlab:help IsolationDistance">IsolationDistance</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help IsolationDistance">IsolationDistance</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'units',
			units = varargin{i+1};
			if ~isimatrix(units,'>=0') | size(units,2) ~= 2,
				error('Incorrect value for property ''units'' (type ''help <a href="matlab:help IsolationDistance">IsolationDistance</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help IsolationDistance">IsolationDistance</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help IsolationDistance">IsolationDistance</a>'' for details).']);
	end
end

% Discard noise and artefacts
artefacts = features(:,3)<=1;
features(artefacts,:) = [];

% By default, process all clusters
if isempty(units),
	units = features(:,[2 3]);
	units = unique(units,'rows');
	units = units(units(:,2)~=0&units(:,2)~=1,:);
end

if strcmp(lower(show),'on'),
	figure;
end
for i = 1:size(units,1),
	unit = units(i,:);
	distance(i,1) = ComputeOnce(features,unit,show);
end

% --------------------------------------------------------------------------------------------------------

function distance = ComputeOnce(features,unit,show)

distance = NaN;
nFeatures = size(features,2)-3;
unitStr = ['(' int2str(unit(1)) ',' int2str(unit(2)) ')'];

% Determine spikes inside vs outside cluster
spikesInCluster = features(:,2)==unit(1) & features(:,3)==unit(2);
nSpikesInCluster = sum(spikesInCluster);
spikesInElectrodeGroup = features(:,2)==unit(1);
nSpikesInElectrodeGroup = sum(spikesInElectrodeGroup);

% We need at least as many spikes as features, but not more than half the total number of spikes
if nSpikesInCluster < nFeatures,
	warning(['Not enough spikes for unit ' unitStr '.']);
	return
end
if nSpikesInCluster > nSpikesInElectrodeGroup/2,
	warning(['Too many spikes for unit ' unitStr '.']);
	return
end

% Compute Mahalanobis distance for spikes inside (mIn) and outside (mOut) the cluster
m = mahal(features(:,4:end),features(spikesInCluster,4:end));
mIn = m(spikesInCluster);
mOut = m(~spikesInCluster);
% Determine the Mahalanobis of the Nth closest spike outside the cluster (where N is the number of spikes inside the cluster)
sOut = sort(mOut);
distance = sOut(nSpikesInCluster);

if strcmp(lower(show),'on'),
	% Compute and plot pdfs
	figure;
	bins = 0:1:300;
	subplot(3,1,1);
	pdfIn = hist(mIn,bins);
	pdfOut = hist(mOut,bins);
	plot(bins,[pdfIn',pdfOut']);
	xlim([0 100]);
	xlabel('Mahalanobis distance (PDF)');
	ylabel('Count');
	legend(['Cluster ' unitStr],'Others','location','southeast');

	% Compute and plot cdfs
	subplot(3,1,2);
	cdfIn = cumsum(pdfIn);
	cdfOut = cumsum(pdfOut);
	plot(bins,[cdfIn',cdfOut']);
	xlim([0 100]);
	xlabel('Mahalanobis distance (CDF)');
	ylabel('Count');
	legend(['Cluster ' unitStr],'Others','location','southeast');

	% Compute and plot ordered distances
	subplot(3,1,3);
	sIn = sort(mIn);
	n = length(sIn);
	plot(1:n,[sIn sOut(1:n)]);
	hold on;
	plot([1 n],[distance distance],'k:');
	xlim([1 n*1.1]);
	xlabel('Spike #');
	ylabel('Mahalanobis distance');
	t = text(n,distance,num2str(round(100*distance)/100));
	set(t,'fontsize',12);
	drawnow;
end
