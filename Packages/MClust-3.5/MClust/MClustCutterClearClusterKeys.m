function MClustCutterClearClusterKeys(figHandle)

% MClustCutterClearClusterKeys(figHandle)

global MClust_Clusters
for iC = 0:length(MClust_Clusters)
    clusterKeys = findobj(figHandle, 'UserData', iC);
    for iK = 1:length(clusterKeys)
        delete(clusterKeys(iK))
    end
end
