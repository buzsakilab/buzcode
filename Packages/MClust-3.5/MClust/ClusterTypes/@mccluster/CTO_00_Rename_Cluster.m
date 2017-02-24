function [MCC, redraw, rekey, undoable] = CTO_RenameCluster(MCC, varargin)

% f = RenameCluster(MCC, varargin)
%
% INPUTS
%     MCC - a MCCluster
% REQUIRES
%
% OUTPUTS
%     MCC - The updated cluster
%
% 
% ncst 26 Nov 02
% ADR 2008
%

extract_varargin;

redraw= true; rekey = true; undoable = true;

newClusterName = inputdlg(['Enter a new name for Cluster ' num2str(iClust) ],...
	'Rename Cluster',1,{MCC.name});
MCC.name = newClusterName{1};
