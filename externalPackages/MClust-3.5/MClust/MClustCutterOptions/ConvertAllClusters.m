function [redraw, rekey, undoable] = ConvertAllClusters

% [redraw, rekey, undoable] = ConvertAllClusters
%
% INPUTS
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
%
% Changes type of iClust
%
% ADR 2003
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
% Extensively modified by ADR to accomodate new ClusterOptions methodology

global MClust_Clusters
global MClust_Directory

availableTypes = dir(fullfile(MClust_Directory, 'ClusterTypes', '@*'));
availableTypes = {availableTypes.name};

msg = {'Select new type.  Clusters will be converted to new type.'};
[Selection, ok] = listdlg('ListString', availableTypes, ...
	'Name', 'Conversion dialog', ...
	'SelectionMode', 'single', ...
	'PromptString', msg);

redraw = false;
rekey = false;
undoable = false;

if ok
	newClass = availableTypes{Selection};
	for iClust = 1:length(MClust_Clusters)
		if ~isa(MClust_Clusters{iClust}, newClass)
			newCluster = feval(newClass(2:end), GetName(MClust_Clusters{iClust}), MClust_Clusters{iClust});
			MClust_Clusters{iClust} = newCluster;
			redraw = true;
			rekey = true;
            undoable = true;
		end
	end
	msgbox('Clusters converted successfully.', 'Convert clusters')
end

