function MCC = mcconvexhull(clusterName, cluster, conversionflag)

% MCC = mcconvexhull(ClusterName, clusterBase)
% Cluster w/convex hulls - no central points
%    
% ADR 2008
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

global MClust_FeatureSources MClust_FeatureNames MCLUST_DEBUG

MCC.xdimNames = {}; 
MCC.ydimNames = {}; 
MCC.xdimSources = {};
MCC.ydimSources = {};
MCC.cx = {}; 
MCC.cy = {};
MCC.name = '';

switch nargin
	case 0
	case 1
		MCC.name = clusterName;
	case 2
        MCC.name = clusterName; % ADR 2008
		if isa(cluster, 'mcconvexhull')
				MCC.xdimNames = cluster.xdimNames;
                MCC.ydimNames = cluster.ydimNames;
                MCC.xdimSources = cluster.ydimSources;
				MCC.ydimSources = cluster.ydimSources;
				MCC.cx = cluster.cx;
				MCC.cy = cluster.cy;
		else
			f = FindInCluster(cluster);
            if ~isempty(f)                
                featuresToUse = 1:length(MClust_FeatureNames);
                timeFeature = strmatch('_Time', MClust_FeatureNames);
                if ~isempty(timeFeature)
                    featuresToUse(timeFeature) = [];
                end
                nDims = length(featuresToUse);
                nBounds = nDims * (nDims-1)/2;
                MCC.xdimNames = cell(nBounds,1);
                MCC.ydimNames = cell(nBounds,1);
                MCC.xdimSources = cell(nBounds,2);
                MCC.ydimSources = cell(nBounds,2);
                iB = 1;
                if false && MCLUST_DEBUG; fig = figure; end                    
                for iC = 1:nDims
                    temp = load(MClust_FeatureSources{iC,1}, '-mat', 'FeatureData');
                    CurrentFeatureDataX = temp.FeatureData(:,MClust_FeatureSources{iC,2});

                    for jC = (iC+1):nDims
                        DisplayProgress(iB, nBounds, 'Title', 'Converting...');
                        temp = load(MClust_FeatureSources{jC,1}, '-mat', 'FeatureData');
                        CurrentFeatureDataY = temp.FeatureData(:,MClust_FeatureSources{jC,2});

                        MCC.xdimNames{iB} = MClust_FeatureNames{featuresToUse(iC)};
                        MCC.xdimSources{iB,1} = MClust_FeatureSources{iC,1};
                        MCC.xdimSources{iB,2} = MClust_FeatureSources{iC,2};

                        MCC.ydimNames{iB} = MClust_FeatureNames{featuresToUse(jC)};
                        MCC.ydimSources{iB,1} = MClust_FeatureSources{jC,1};
                        MCC.ydimSources{iB,2} = MClust_FeatureSources{jC,2};

                        k = convexhull(CurrentFeatureDataX(f), CurrentFeatureDataY(f));
                        MCC.cx{iB} = CurrentFeatureDataX(f(k));
                        MCC.cy{iB} = CurrentFeatureDataY(f(k));
                        
                        if false && MCLUST_DEBUG
                            figure(fig); clf;
                            plot(CurrentFeatureDataX(f), CurrentFeatureDataY(f), '.');
                            hold on
                            plot(MCC.cx{iB}, MCC.cy{iB}, 'k');
                            drawnow; 
                        end                                                     
                        
                        iB = iB + 1;                        

                    end
                end
            end
			DisplayProgress close;
        end
    case 3 % convert from another dataset
        MCC.name = clusterName; % ADR 2008
        if isa(cluster, 'mcconvexhull')
            nBounds = length(cluster.xdimNames);
            MCC.xdimNames = cluster.xdimNames;
            MCC.cx = cluster.cx;
            MCC.ydimNames = cluster.ydimNames;
            MCC.cy = cluster.cy;

            MCC.xdimSources = cell(nBounds,2);
            MCC.ydimSources = cell(nBounds,2);

            for iB = 1:nBounds
                featureX = strmatch(cluster.xdimNames{iB}, MClust_FeatureNames);
                if (length(featureX) ~= 1)
                    errordlg(['Cannot find feature match for ' cluster.xdimNames{iB}], ...
                        'MClust error', 'modal');
                    MCC = [];
                    return;
                end
                MCC.xdimSources{iB,1} = MClust_FeatureSources{featureX,1};
                MCC.xdimSources{iB,2} = MClust_FeatureSources{featureX,2};
                
                featureY = strmatch(MCC.ydimNames{iB}, MClust_FeatureNames);
                if (length(featureY) ~= 1)
                    errordlg(['Cannot find feature match for ' cluster.ydimNames{iB}], ...
                        'MClust error', 'modal');
                    MCC = [];
                    return;
                end
                MCC.ydimSources{iB,1} = MClust_FeatureSources{featureY,1};
                MCC.ydimSources{iB,2} = MClust_FeatureSources{featureY,2};
            end
        else 
            MCC = [];
        end

	otherwise
		error('Incorrect number of inputs to MCConvexHull');
end

MCC = class(MCC, 'mcconvexhull');