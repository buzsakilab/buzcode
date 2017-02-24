function MClustCutterRedrawClusterKeys(figHandle, startCluster)

% MClustCutterRedrawClusterKeys(figHandle, startCluster)

global MClust_Clusters MClust_Colors MClust_Directory

if nargin == 1
    sliderHandle = findobj(figHandle, 'Tag', 'ScrollClusters');
    startCluster = floor(-get(sliderHandle, 'Value'));
end

endCluster = floor(min(startCluster + 16, length(MClust_Clusters)));

ClusterOptionsFiles = FindFiles('*.m','CheckSubdirs',0, ...
	'StartingDirectory', fullfile(MClust_Directory, 'ClusterOptions'));
Extra_Options = cell(length(ClusterOptionsFiles),1);
for iCOF = 1:length(ClusterOptionsFiles)
	[dummy_fd Extra_Options{iCOF} ext] = fileparts(ClusterOptionsFiles{iCOF});
end

for iC = startCluster:endCluster
	
	if iC == 0
		CreateClusterKeys(figHandle, iC, 0.35, 0.9 - 0.05 * (iC - startCluster), 'MClustCutterCallbacks');
	else % iC > 0
		ClusterOptionsFiles = FindFiles('CTO_*.m','CheckSubdirs',0, ...
			'StartingDirectory', fullfile(MClust_Directory, 'ClusterTypes', ['@' class(MClust_Clusters{iC})]));
		ClusterType_Options = cell(length(ClusterOptionsFiles),1);
		for iCOF = 1:length(ClusterOptionsFiles)
			[dummy_fd ClusterType_Options{iCOF} ext] = fileparts(ClusterOptionsFiles{iCOF});
			ClusterType_Options{iCOF} = ClusterType_Options{iCOF}(5:end);
		end
		
		CreateClusterKeys(figHandle, iC, 0.35, 0.9 - 0.05 * (iC - startCluster), 'MClustCutterCallbacks', ...
			ClusterType_Options{:}, '--------------', Extra_Options{:});

		if iC+1 > size(MClust_Colors,1)
			MClust_Colors(end+1,:) = MClust_Colors(end,:);
		end
	end
end
