function MClustCutterRedrawAxes(figHandle, varargin)

global MClust_Clusters MClust_Colors MClust_Hide MClust_UnaccountedForOnly 
global MClust_ClusterIndex MClust_FeatureData 

global MClust_CurrentFeatures % used to keep track of which features are currently in memory
global MClust_CurrentFeatureNames % 
global MClust_CurrentFeatureData
global MClust_FeatureNames % names of features
global MClust_FeatureSources % <filenames, number pairs> for finding features in fd files

global MClust_ClusterCutWindow_Marker
global MClust_ClusterCutWindow_MarkerSize
global MClust_CHDrawingAxisWindow_Pos;

% -- get variables
full = 0;
extract_varargin;

nClust = length(MClust_Clusters);

drawingFigHandle = findobj('Type', 'figure', 'Tag', 'CHDrawingAxisWindow');  % figure to draw in

xdimHandle = findobj(figHandle, 'Tag', 'xdim');
xdim = get(xdimHandle, 'Value');           % x dimemsion to plot
ydimHandle = findobj(figHandle, 'Tag', 'ydim');  
ydim = get(ydimHandle, 'Value');           % y dimension to plot
markerHandle = findobj(figHandle, 'Tag', 'PlotMarker');
markerString = get(markerHandle, 'String');
markerValue = get(markerHandle, 'Value');
MClust_ClusterCutWindow_Marker = markerValue;
marker = markerString{markerValue};
markerSizeHandle = findobj(figHandle, 'Tag', 'PlotMarkerSize');
markerSizeString = get(markerSizeHandle, 'String');
markerSizeValue = get(markerSizeHandle, 'Value');
MClust_ClusterCutWindow_MarkerSize = markerSizeValue;
markerSize = str2double(markerSizeString{markerSizeValue});

% converted back to work by disk access (ADR 2008 - turns out this is faster
% get xdim
if (MClust_CurrentFeatures(1) ~= xdim)
    temp = load(MClust_FeatureSources{xdim,1}, '-mat', 'FeatureData');
    MClust_CurrentFeatureData(:,1) = temp.FeatureData(:,MClust_FeatureSources{xdim,2});
    MClust_CurrentFeatures(1) = xdim;
    MClust_CurrentFeatureNames{1} = MClust_FeatureNames{xdim};
end
% get ydim
if (MClust_CurrentFeatures(2) ~= ydim)
    temp = load(MClust_FeatureSources{ydim,1}, '-mat', 'FeatureData');
    MClust_CurrentFeatureData(:,2) = temp.FeatureData(:,MClust_FeatureSources{ydim,2});
    MClust_CurrentFeatures(2) = ydim;
    MClust_CurrentFeatureNames{2} = MClust_FeatureNames{ydim};
end

if isempty(drawingFigHandle)
    % create new drawing figure
    drawingFigHandle = ...
        figure('Name', 'Cluster Cutting Window',...
        'NumberTitle', 'off', ...
        'Tag', 'CHDrawingAxisWindow', ...
        'KeyPressFcn', 'MClustCutterKeyPress','Position',MClust_CHDrawingAxisWindow_Pos);
else
    % figure already exists -- select it
    figure(drawingFigHandle);
end

% have to a complete redraw
if ~full
    curAxis = axis;
end
clf;
hold on;
if full %%% Added by JCJ Aug 2007 to stabilize redraw of display
    set(gca, 'XLim', [min(MClust_CurrentFeatureData(:,1)) max(MClust_CurrentFeatureData(:,1))+0.0001]);
    set(gca, 'YLim', [min(MClust_CurrentFeatureData(:,2)) max(MClust_CurrentFeatureData(:,2))+0.0001]);	
else
    axis(curAxis);
end
for iC = 0:nClust
    if ~MClust_Hide(iC+1)
        HideClusterHandle = findobj(figHandle, 'UserData', iC, 'Tag', 'HideCluster');
        if iC == 0
            if MClust_UnaccountedForOnly
                MClust_ClusterIndex = ProcessClusters(MClust_CurrentFeatureData, MClust_Clusters);
                f = (MClust_ClusterIndex == 0);
                 figure(drawingFigHandle);
                h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
            else
                figure(drawingFigHandle);
                h = plot(MClust_CurrentFeatureData(:,1), MClust_CurrentFeatureData(:,2), marker);
            end
        else         
            [f,MClust_Clusters{iC}] = FindInCluster(MClust_Clusters{iC});
            if isempty(f) && ~isempty(HideClusterHandle)
                set(HideClusterHandle, 'Enable', 'off');
            else 
                set(HideClusterHandle, 'Enable', 'on');
            end
            figure(drawingFigHandle);
            h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
        end
        set(h, 'Color', MClust_Colors(iC+1,:));
        set(h, 'Tag', 'ClusterLine', 'UserData', iC);
        set(h, 'MarkerSize', markerSize);
		if iC > 0
			try
				DrawOnAxis(MClust_Clusters{iC}, xdim, ydim, MClust_Colors(iC+1,:), gca); 
			end
		end
    end
end
figure(drawingFigHandle);
if full
    set(gca, 'XLim', [min(MClust_CurrentFeatureData(:,1)) max(MClust_CurrentFeatureData(:, 1))+0.0001]);
    set(gca, 'YLim', [min(MClust_CurrentFeatureData(:,2)) max(MClust_CurrentFeatureData(:, 2))+0.0001]);
else
    axis(curAxis);
end
xlabel(MClust_CurrentFeatureNames{1},'interpreter','none');
ylabel(MClust_CurrentFeatureNames{2},'interpreter','none');
zoom on

contourWindow = findobj('Type', 'figure', 'Tag', 'ContourWindow');
if ~isempty(contourWindow)
    mkContours(drawingFigHandle, 'figHandle', contourWindow);
end
figure(drawingFigHandle);
