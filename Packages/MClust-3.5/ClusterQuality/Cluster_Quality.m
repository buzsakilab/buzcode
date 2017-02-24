function [CluSep, m] = Cluster_Quality(Fet, ClusterSpikes)

% [CluSep, m] = Cluster_Quality(Fet, ClusterSpikes, doplot)
%
% Calculate measures of cluster quality: L-ratio and Isolation Distance
%
% Inputs:   Fet:           N by D array of feature vectors (N spikes, D dimensional feature space)
%           ClusterSpikes: Index into Fet which lists spikes from the cell whose quality is to be evaluated.

% calculate square mahalanobis distance from the center of the cluster to
% all other spikes

if length(ClusterSpikes) < size(Fet,2)
    disp(['  Insufficient spikes to calculate Mahalanobis distance: at least ' num2str(size(Fet,2)) ' spikes required.'])
    
    ClusSep.nSpikes = length(ClusterSpikes);
    CluSep.IsolationDist = NaN;
    CluSep.L  = NaN;
    CluSep.Lratio = NaN;
    CluSep.df = NaN;
    return
end

m = mahal(Fet,Fet(ClusterSpikes,:));

% number of spikes in the cluster
ClusSep.nSpikes = length(ClusterSpikes);

% Get Isolation Distance
CluSep.IsolationDist = IsolationDistance(Fet,ClusterSpikes,m);

% Get L_ratio
[CluSep.L CluSep.Lratio CluSep.df] = L_ratio(Fet,ClusterSpikes,m);