function MCC = mccluster(clusterName, cluster)

% MCC = mccluster(ClusterName, clusterBase)
% Cluster
%    
% ADR 2008
% PL 2000-2002
% cowen 2002 added a restriction -- ForbiddenPoints - great if you need to exclude points cut out by non hull bound criteria.
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

global MClust_FeatureSources

MCC.name = '';
MCC.xdimNames = {}; 
MCC.ydimNames = {}; 
MCC.xdimSources = {};
MCC.ydimSources = {};
MCC.cx = {}; 
MCC.cy = {};
MCC.AddFlag = [];
MCC.recalc = -1;
MCC.myPoints = [];
MCC.myOrigPoints = [];
MCC.ForbiddenPoints = []; % cowen 2002 restrict results to a subset, regardless of convex hulls.


switch nargin
	case 0		
	case 1
		if isa(clusterName, 'mccluster')
			MCC = clusterName; % loading from copy
		else
			MCC.name = clusterName;
		end
	case 2
		MCC.myPoints = FindInCluster(cluster);
		MCC.myOrigPoints = FindInCluster(cluster);
		MCC.name = clusterName;
	otherwise
		error('Incorrect number of inputs to MCConvexHull');
end

MCC.recalc = 1;

MCC = class(MCC, 'mccluster');