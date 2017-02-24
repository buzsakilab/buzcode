function [redraw, rekey, undoable] = DeleteCluster(iClust)

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
global MClust_Colors

redraw= true; rekey = true; undoable = true;
    
figHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');

MClustCutterClearClusterKeys(figHandle);
MClust_Clusters(iClust) = [];
MClust_Hide(iClust+1) = [];
MClust_Colors(iClust+1,:) = [];
	
warndlg(['Cluster ' num2str(iClust) ' deleted.']);
