function [redraw, rekey, undoable] = CopyCluster(iClust)

% f = CopyCluster(iClust, varargin)
%
% INPUTS
%     iClust - index into MClust_Clusters
% REQUIRES
%  figHandle from parent
%
% OUTPUTS
%     MCC - The updated cluster
%
% 
% ncst 26 Nov 02
% ADR 2008
%

global MClust_Clusters
global MClust_Hide

redraw= true; rekey = true; undoable = true;

newCluster = MClust_Clusters{iClust};
newCluster = SetName(newCluster, ['Copy of ' GetName(newCluster)]);
MClust_Clusters{end+1} = newCluster;
MClust_Hide(end+1) = 0;
	
warndlg(['Cluster ' num2str(iClust) ' copied to cluster ' ...
	num2str(length(MClust_Clusters)) '.'], 'Copy successful');
