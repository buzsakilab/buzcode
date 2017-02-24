function [redraw, rekey, undoable] = ConvertCluster(iClust)

% [redraw, rekey, undoable] = Convertcluster(iClust)
%
% INPUTS
%    iClust
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

oldType = class(MClust_Clusters{iClust});

availableTypes = dir(fullfile(MClust_Directory, 'ClusterTypes', '@*'));
availableTypes = {availableTypes.name};

msg = {sprintf('%s: (%s).', GetName(MClust_Clusters{iClust}), class(MClust_Clusters{iClust})),...
	'Select new type.  Cluster will be converted to new type.'};
[Selection, ok] = listdlg('ListString', availableTypes, ...
	'InitialValue', strmatch(class(MClust_Clusters{iClust}), availableTypes), ...
	'Name', 'Conversion dialog', ...
	'SelectionMode', 'single', ...
	'PromptString', msg);

redraw = false; rekey = false; undoable = false;

if ok
	newClass = availableTypes{Selection};
	if ~isa(MClust_Clusters{iClust}, newClass)
		newCluster = feval(newClass(2:end), GetName(MClust_Clusters{iClust}), MClust_Clusters{iClust});
		MClust_Clusters{iClust} = newCluster;
		redraw = true;
		rekey = true;
        undoable = true;
	end
end

