function [redraw, rekey, undoable] = RunWaveformCutter(iClust)
% RunWaveformCutter(iClust)
%
% INPUTS
%    iClust
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder

% ADR 2003
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.
% Extensively modified by ADR to accomodate new ClusterOptions methodology
% 2008 Mar - ADR undo is done in cutter.  Not here.

redraw = false; rekey = false; undoable = false; 

global MClust_Clusters MClust_TTData MClust_ChannelValidity

% check to make sure cluster is compatible
if isempty(which(fullfile(class(MClust_Clusters{iClust}), 'Restrict_Points')))
	errordlg(sprintf('Cluster [%s] is of type [%s], which does not support Restrict_Points, necessary for WaveformCutter.', ...
		GetName(MClust_Clusters{iClust}), class(MClust_Clusters{iClust})));
	return;
end

subf = FindInCluster(MClust_Clusters{iClust});

wv = ExtractCluster(MClust_TTData, subf);

t  = Range(wv,'ts');
wv = Data(wv);
WaveformCutter(t,wv,MClust_ChannelValidity,iClust,subf);

