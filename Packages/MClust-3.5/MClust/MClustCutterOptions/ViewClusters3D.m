function [redraw, rekey, undoable] = ViewClusters3D

% [redraw, rekey, undoable] = ViewClusters3D
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
% ADR 2008
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

redraw = false; rekey = false; undoable = false;

global MClust_Clusters MClust_Hide MClust_Colors MClust_UnaccountedForOnly
global MClust_FeatureNames MClust_CurrentFeatureData MClust_FeatureSources

figHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');
xdimHandle = findobj(figHandle, 'Tag', 'xdim');
xdim = get(xdimHandle, 'Value');           % x dimemsion to plot
ydimHandle = findobj(figHandle, 'Tag', 'ydim');
ydim = get(ydimHandle, 'Value');           % y dimension to plot
for iC = 1:length(MClust_FeatureNames)
	if ~(iC == xdim) && ~(iC==ydim)
		break
	end
end
[Selection, OK] = listdlg('ListString', MClust_FeatureNames, 'SelectionMode', 'single', 'Name', 'zdim',...
	'InitialValue', iC, ...
	'PromptString', {sprintf('Using (%s x %s).',MClust_FeatureNames{xdim}, MClust_FeatureNames{ydim}),...
	'Select z dimension.'});
zdim = Selection;
if ~OK
	return
end

markerHandle = findobj(figHandle, 'Tag', 'PlotMarker');
markerString = get(markerHandle, 'String');
markerValue = get(markerHandle, 'Value');
marker = markerString{markerValue};
markerSizeHandle = findobj(figHandle, 'Tag', 'PlotMarkerSize');
markerSizeString = get(markerSizeHandle, 'String');
markerSizeValue = get(markerSizeHandle, 'Value');
markerSize = str2double(markerSizeString{markerSizeValue});


temp = load(MClust_FeatureSources{zdim,1}, '-mat', 'FeatureData');
CurrentFeatureDataZ = temp.FeatureData(:,MClust_FeatureSources{zdim,2});

drawingFigHandle = figure('Name', 'Clusters 3D',...
	'NumberTitle', 'off', ...
	'Tag', 'DrawingAxisWindow3D');
clf; hold on;
nClust = length(MClust_Clusters);
for iC = 0:nClust
	if ~MClust_Hide(iC+1)
		if iC == 0
			if MClust_UnaccountedForOnly
				MClust_ClusterIndex = ProcessClusters(MClust_CurrentFeatureData, MClust_Clusters);
				f = (MClust_ClusterIndex == 0);
				figure(drawingFigHandle);
				h = plot3(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), CurrentFeatureDataZ(f),'.',...
					'marker', marker, 'markerSize', markerSize, 'color', MClust_Colors(iC+1,:));
			else
				figure(drawingFigHandle);
				h = plot3(MClust_CurrentFeatureData(:,1), MClust_CurrentFeatureData(:,2), CurrentFeatureDataZ(:),'.',...
					'marker', marker, 'markerSize', markerSize, 'color', MClust_Colors(iC+1,:));
			end
		else
			[f,MClust_Clusters{iC}] = FindInCluster(MClust_Clusters{iC});
			figure(drawingFigHandle);
			h = plot3(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), CurrentFeatureDataZ(f),'.',...
				'marker', marker, 'markerSize', markerSize, 'color', MClust_Colors(iC+1,:));
		end
	end
end

xlabel(MClust_FeatureNames{xdim});
ylabel(MClust_FeatureNames{ydim});
zlabel(MClust_FeatureNames{zdim});
rotate3d on
view(3);


