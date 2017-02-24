function [exploration,sws,rem] = BrainStates(s,t,f,q,emg,varargin)

%BrainStates - Determine brain state using LFP, EMG and movement.
%
%  USAGE
%
%    [exploration,sws,rem] = BrainStates(s,t,f,q,emg,<options>)
%
%    s              spectrogram
%    t              time bins for spectrogram
%    f              frequency bins for spectrogram
%    q              'quiescence' vector, obtained using <a href="matlab:help QuietPeriods">QuietPeriods</a> (see NOTE)
%    emg            electromyogram (pass [] if missing)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nClusters'   number of clusters for K-means clustering (default = 2)
%     'method'      whether brain rhythms should be determined based on a PCA
%                   on the spectrogram ('pca'), or using predefined frequency
%                   band ratios, such as the theta/delta ratio ('direct', by
%                   default) or the heuristic ratios described in Gervasoni et
%                   al. 2004 ('ratios')
%     'nComponents' number of principal components (for spectrogram PCA)
%                   (default = automatically computed to account for 85% of
%                   the variance)
%     'show'        plot K-means clustering in feature space ('kmeans'),
%                   clusters on spectrogram ('clusters') or both ('all')
%                   (default = 'none')
%    =========================================================================
%
%  OUTPUT
%
%    exploration    at each time t, whether the rat was exploring
%    sws            at each time t, whether the rat was in slow wave sleep
%    rem            at each time t, whether the rat was in rapid eye movement sleep
%
%  NOTE
%
%    Quiescence q is a time series of boolean values specifying for each
%    timestamp whether the animal was quiet (as determined by <a href="matlab:help QuietPeriods">QuietPeriods</a>).
%
%  SEE
%
%    See also QuietPeriods.

% Copyright (C) 2008-2011 by Micha??l Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
nClusters = 2;
window = 5*1250;
show = 'none';
nComponents = 0;
method = 'direct';

% Check number of parameters
if nargin < 5,
	error('Incorrect number of parameters (type ''help <a href="matlab:help BrainStates">BrainStates</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help BrainStates">BrainStates</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help BrainStates">BrainStates</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'nclusters',
			nClusters = varargin{i+1};
			if ~isiscalar(nClusters,'>0'),
				error('Incorrect value for property ''nClusters'' (type ''help <a href="matlab:help BrainStates">BrainStates</a>'' for details).');
			end
		case 'method',
			method = lower(varargin{i+1});
			if ~isstring_FMAT(method,'pca','direct','ratios'),
				error('Incorrect value for property ''method'' (type ''help <a href="matlab:help BrainStates">BrainStates</a>'' for details).');
			end
		case 'ncomponents',
			nComponents = varargin{i+1};
			if ~isiscalar(nComponents,'>0'),
				error('Incorrect value for property ''nComponents'' (type ''help <a href="matlab:help BrainStates">BrainStates</a>'' for details).');
			end
		case 'show',
			show = lower(varargin{i+1});
			if ~isstring_FMAT(show,'kmeans','clusters','all','none'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help BrainStates">BrainStates</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help BrainStates">BrainStates</a>'' for details).']);
	end
end

% Compute and interpolate EMG and quiescence at spectrogram timestamps
emg0 = [];
if ~isempty(emg),
	emg0 = Interpolate(emg,t,'trim','off');
	emg0 = Smooth(abs(emg0(:,2)),5);
end
q0 = [];
if all(q(:,2)==0|q(:,2)==1),
	q0 = Interpolate(double(q),t,'trim','off');
	q0 = round(q0(:,2));
else
	q0 = Interpolate(q,t,'trim','off');
	q0 = q0(:,2);
end

% No movement info => undetermined
usable = ~isnan(q0);

% Exploration corresponds to movement periods
exploration = zeros(size(t));
exploration(usable) = ~q0(usable);

% Sleep/rest corresponds to quiescence
sleep = zeros(size(t));
sleep(usable) = q0(usable);
sleep = logical(sleep);

% Get theta/delta ratio
bands = SpectrogramBands(s,f);

% Determine features for automatic data clustering
switch(method),
	case 'pca',
		% Compute PCA on spectrogram and reduce dimensionality
		S = s(:,sleep)';
		S(:,f>30) = 0;
		[eigenvectors,projected,lambda] = princomp(S,'econ');
		if nComponents == 0,
			nComponents = find(cumsum(lambda)/sum(lambda)>0.85);
			nComponents = nComponents(1);
		end
		eigenvectors = eigenvectors(:,1:nComponents);
		lambda = lambda(1:nComponents);
		projected = projected(:,1:nComponents);
		features = projected;
	case 'direct',
		% Use theta/delta ratio
		features = bands.ratio(sleep);
	case 'ratios',
		% Use heuristic ratios
		features = [bands.ratio1(sleep) bands.ratio2(sleep)];
end

% Add EMG to feature list
if ~isempty(emg0),
	features = [features emg0(sleep)];
end

% Normalize columns (we need this because K-means clustering uses Euclidean
% distances)
features = zscore(features);

% Cluster (K-means)
cluster = kmeans(features,nClusters);

% Identify clusters
sws = logical(zeros(size(t)));
rem = logical(zeros(size(t)));
% 1) Measure mean theta/delta ratio for each state
r = bands.ratio(sleep);
for i = 1:nClusters,
	ratio(i) = mean(r(cluster==i));
end
% 2) REM is the quiet state with highest theta/delta ratio
[unused,i] = sort(ratio(:),1,'descend');
h = i(1);
rem(sleep) = cluster==h;
highest = ratio(h);
% 3) SWS is the quiet state with lowest theta/delta ratio
[unused,i] = sort(ratio(:),1,'ascend');
l = i(1);
sws(sleep) = cluster==l;
lowest = ratio(l);
% 4) REM must have a theta/delta ratio at least twice as high as SWS
if highest < 2*lowest,
	sws(sleep) = 1;
	rem(sleep) = 0;
end

% Show K-means clustering in feature space
if strcmp(show,'kmeans') | strcmp(show,'all'),
	figure;
	% Determine number N of subplots (depends on number of principal components kept; clip at 15)
	n = size(features,2);
	N = min([n*(n-1)/2 21]);
	if N == 0, N = 1; end
	m0 = round(sqrt(N));
	n0 = ceil(N/m0);
	colors = jet(nComponents);
	if strcmp(method,'pca'),
		% Plot principal components
		subplot(m0+1,1,1);hold on;
		l = 'h = legend(';
		for i = 1:nComponents,
			plot(f,eigenvectors(:,i),'color',colors(i,:));
			l = [l '''' int2str(i) ''','];
		end
		l = [l(1:end-1) ',''location'',''northoutside'',''orientation'',''horizontal'');'];
		eval(l);
		set(h,'color',[1 1 1]*.7);
		shift = 1;
	else
		shift = 0;
	end
	% List pairs of principal component IDs to display in successive 2D projections
	colors = jet(nClusters);
	if N > 1,
		feature1 = [];
		feature2 = [];
		i = 1;
		while length(feature1) < N,
			feature1 = [feature1 repmat(i,1,n-i)];
			feature2 = [feature2 i+1:n];
			i = i+1;
		end
		% Plot 2D projections of data in PC space, colored by cluster ID
		for i = 1:N,
			subplot(m0+shift,n0,i+shift*n0);hold on;
			for j = 1:size(features,1),
				plot(features(j,feature1(i)),features(j,feature2(i)),'.','color',colors(cluster(j),:));
			end
			xlabel([int2str(feature1(i)) ' vs ' int2str(feature2(i))]);
		end
	else
		subplot(shift+1,1,shift+1);hold on;
		r = rand(size(features));
		for j = 1:size(features,1),
			plot(r(j),features(j),'.','color',colors(cluster(j),:));
		end
	end
end

% Show clusters on spectrogram
if strcmp(show,'clusters') | strcmp(show,'all'),
	figure;
	% 1) Plot spectrogram and theta/delta ratio
	subplot(2,1,1);hold on;
	PlotColorMap(log(s),1,'cutoffs',[0 12],'x',t,'y',f);
	plot(t,bands.ratio,'k');
	l = 'h = legend(''\theta/\delta ratio'',''q''';
	ylim([0 30]);
	yLim = ylim;
	dy = yLim(2) - yLim(1);
	y0 = yLim(1);
	% ... quiescence...
	plot(t,5*ZeroToOne(q0)+y0+0.5*dy,'b');
	% ... and EMG
	if ~isempty(emg0),
		plot(t,ZeroToOne(emg0)*5+y0+0.65*dy,'r');
		l = [l ',''emg'''];
	end
	l = [l ');'];
	eval(l);
	set(h,'color',[1 1 1]);
	xlim([t(1) t(end)]);
	% 2) Plot spectrogram and clusters
	subplot(2,1,2);hold on;
	PlotColorMap(log(s),1,'cutoffs',[0 12],'x',t,'y',f);
	ylim([0 30]);
	l = 'h = legend(';
	yLim = ylim;
	dy = yLim(2) - yLim(1);
	y0 = yLim(1);
	% Plot clustering output
	c = nan(size(t));
	c(sleep) = cluster-1;
	plot(t,5*ZeroToOne(c)+y0+0.35*dy,'k');
	l = [l '''cluster'',''location'',''north'',''orientation'',''horizontal'');'];
	eval(l);
	% Show SWS/REM
	if ~isempty(sws),
		z = nan(size(t));
		z(sws) = t(sws);
		plot(z,15*~isnan(z),'b','linewidth',5);
	end
	if ~isempty(rem),
		z = nan(size(t));
		z(rem) = t(rem);
		plot(z,15*~isnan(z),'r','linewidth',5);
	end
	xlim([t(1) t(end)]);
	set(h,'color',[1 1 1]);
end
