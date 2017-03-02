function IsolDist = IsolationDistance(Fet, ClusterSpikes, m)

% IsolDist = IsolationDistance(Fet, ClusterSpikes, m)
%
% Isolation Distance
% Measure of cluster quality
%add
% Inputs:   Fet:           N by D array of feature vectors (N spikes, D dimensional feature space)
%           ClusterSpikes: Index into Fet which lists spikes from the cell whose quality is to be evaluated.
%           m:             squared mahalanobis distances, default is to
%                          calculate them directly
%
% Created by Ken Harris

% find # of spikes in this cluster
if nargin < 3
	nSpikes = size(Fet,1);
else
	nSpikes = size(m,1);
end

nClusterSpikes = length(ClusterSpikes);
InClu = ismember(1:nSpikes, ClusterSpikes);

% check there are enough spikes (but not more than half)
if length(ClusterSpikes) < size(Fet,2) | length(ClusterSpikes)>nSpikes/2
	IsolDist = NaN;
	return
end

% mark spikes which are not cluster members
NoiseSpikes = setdiff(1:nSpikes, ClusterSpikes);

%%%%%%%%%%% compute mahalanobis distances %%%%%%%%%%%%%%%%%%%%%
if nargin < 3
	m = mahal(Fet, Fet(ClusterSpikes,:));
end

mCluster = m(ClusterSpikes); % mahal dist of spikes in the cluster
mNoise = m(NoiseSpikes); % mahal dist of all other spikes

% calculate point where mD of other spikes = n of this cell
if (nClusterSpikes < nSpikes/2)
	[sorted order] = sort(mNoise);
	IsolDist = sorted(nClusterSpikes);
else
	IsolDist = NaN; % If there are more of this cell than every thing else, forget it.
end