function [redraw, rekey, undoable] = CheckCluster(iClust)
% [redraw, rekey, undoable] = CheckCluster(iClust)
%
% INPUTS
%    iClust
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
%
% ADR 2008
%
% Status: PROMOTED (Release version)
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.
% Extensively modified by ADR to accomodate new ClusterOptions methodology

redraw = false; rekey = false; undoable = false; % don't need to update

global MClust_TTData MClust_Clusters MClust_FeatureSources
global MClust_TTdn MClust_TTfn MClust_TText
global MClust_FeatureTimestamps MClust_ChannelValidity

run_checkclust = 1;

[f MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});
[L_Extra,L_Ratio,IsolationDist,Dists] = ClusterSeparation(f,...
	fullfile(MClust_TTdn, [MClust_TTfn MClust_TText]),MClust_ChannelValidity,iClust);

if isempty(f) 
	run_checkclust = 0;
	msgbox('No points in cluster.')
	return
end

save_jpegs = 0;
if iClust == 0
	clustTT = MClust_TTData;
else
	[clustTT, was_scaled] = ExtractCluster(MClust_TTData, f);
end
[curdn,curfn] = fileparts(MClust_TTfn);
if was_scaled
	title_string = [ 'Only ' num2str(length(Range(clustTT,'ts'))) ' of ' num2str(was_scaled) ' spikes shown.'];
else
	title_string = [];
end
ClustInd = repmat(0,size(MClust_FeatureTimestamps));
ClustInd(f) = 1;
CheckCluster([curfn, '--Cluster', num2str(iClust)], clustTT, MClust_FeatureTimestamps(1), MClust_FeatureTimestamps(end), save_jpegs,...
	title_string, was_scaled,'L_Extra',L_Extra,'L_Ratio',L_Ratio,'IsolationDist',IsolationDist,...
	'Dists',Dists,'ClustInd',ClustInd);


