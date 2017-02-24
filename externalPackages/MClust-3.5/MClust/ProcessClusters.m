function clusterIndex = ProcessClusters(data, clusters)

% clusterIndex = ProcessClusters(data, clusters)
%
% INPUT:
%   data - MClust feature data
%   clusters - MClust cluster objects
% OUTPUT:
%   index showing to which cluster each spike belongs
%   spikes belonging to multiple clusters get assigned the error signal -1
%
% ADR 1999
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


[nSamps, nDims] = size(data);
clusterIndex = zeros(nSamps,1);
for iC = 1:length(clusters)
   f = FindInCluster(clusters{iC});
   f2 = (clusterIndex(f)>0);
   clusterIndex(f) = iC;
   clusterIndex(f2) = -1;
end


