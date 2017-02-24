function KlustaKwikCallbacks(varargin)

% ViewClustersCallbacks
%
% Callbacks for view clusters window
%
% ADR 1998
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

%JCJ Feb 2008 -- Added intermediate step to KKwik export to eliminate
%                function indexing error (MATLAB 7.3.0)


%------------------------------------

% globals 
global MClust_Clusters MClust_FeatureSources MClust_FeaturesToUse MClust_FeatureTimestamps
global MClust_TTData MClust_TTfn MClust_ChannelValidity

global KlustaKwik_Clusters
global KKClust KlustaKwik_AllHold

global MClust_KKDecisionWindow_Pos
global MClust_KK2D_Pos 
global MClust_KK3D_Pos
global MClust_KKContour_Pos

global MClust_Colors MClust_FeatureNames
      
% get the various handles
KKDW = findobj('Tag', 'KKDecisionWindow');
if ~isempty(KKDW)
    MClust_KKDecisionWindow_Pos = get(KKDW,'Position');
end

drawingFigHandle2D = findobj('Type', 'figure', 'Tag', 'KK2D');
if ~isempty(drawingFigHandle2D)
    MClust_KK2D_Pos = get(drawingFigHandle2D,'Position');
end

drawingFigHandle3D = findobj('Type', 'figure', 'Tag', 'KK3D');
if ~isempty(drawingFigHandle3D)
    MClust_KK3D_Pos = get(drawingFigHandle3D,'Position');
end

drawingFigHandleContour = findobj('Tag', 'KKContour');
if ~isempty(drawingFigHandleContour)
    MClust_KKContour_Pos = get(drawingFigHandleContour,'Position');
end

if ~isempty(varargin)
  figHandle = findobj('Tag', 'KKDecisionWindow');
  callbackTag = varargin{1};
  if length(varargin) > 1
      cboHandle = varargin{2};
  else
      cboHandle = [];
  end
else  
  cboHandle = gcbo;
  figHandle = gcf;
  callbackTag = get(cboHandle, 'Tag');
end

%--------------------------------------
% main switch
switch callbackTag
    
%%
	case 'ExportClustersSeparate'
		
	clusterStart = length(MClust_Clusters)+1;
    for iC = 1:length(KlustaKwik_Clusters)
        KeepYN = get(findobj(figHandle, 'Tag', 'KeepCluster', 'UserData', iC), 'Value');
        Color = get(findobj(figHandle, 'Tag', 'ChooseColor', 'UserData', iC), 'BackgroundColor');
        if KeepYN
            MClust_Clusters{end+1} = precut(['KKwik ' num2str(iC)]);
			MClust_Colors(1+length(MClust_Clusters),:) = Color;
			MClust_Clusters{end} = AddIndices(MClust_Clusters{end}, ...
				FindInCluster(KlustaKwik_Clusters{iC}));
        end%if   
    end%for
	
	% get type to export as
	exportTypeHandle = findobj('Tag','ExportAsType');
	exportType = get(exportTypeHandle, 'String');
	exportType = exportType{get(exportTypeHandle, 'Value')};
	for iC = clusterStart:length(MClust_Clusters)
		if ~isa(MClust_Clusters{iC}, exportType);
			MClust_Clusters{iC} = feval(exportType(2:end), GetName(MClust_Clusters{iC}), MClust_Clusters{iC});
		end
	end
  
	% close KKwik windows
	% CloseKK() - this is a pain, don't close ADR Feb 2008

   % start manual cut or update if already open
   CutterFig = findobj('Type','figure','Tag', 'ClusterCutWindow');
   if isempty(CutterFig)
	   MClustCutter();
   else
	   MClustCutterClearClusterKeys(CutterFig);
	   MClustCutterRedrawClusterKeys(CutterFig);
   end

case 'ExportClustersMerge'

	clusterStart = length(MClust_Clusters)+1; % ADR 2008
    nC = length(KlustaKwik_Clusters);
	done = zeros(nC,1);	
	keepYN = zeros(nC,1);
	color = zeros(nC,3);
	clusterTo = zeros(nC,1);
	for iC = 1:nC;
		keepYN(iC) = get(findobj(figHandle, 'Tag', 'KeepCluster', 'UserData', iC), 'Value');
		color(iC,:) = get(findobj(figHandle, 'Tag', 'ChooseColor', 'UserData', iC), 'BackgroundColor');
	end
	
	for iC = 1:nC
		if ~done(iC) && keepYN(iC)
			MClust_Clusters{end+1} = precut(['KKwik ' num2str(iC)]);
			MClust_Colors(1+length(MClust_Clusters),:) = color(iC,:);
			for jC = iC:nC
				if all(color(iC,:)==color(jC,:))
					MClust_Clusters{end} = AddIndices(MClust_Clusters{end}, ...
						FindInCluster(KlustaKwik_Clusters{jC})); 
					done(jC) = 1;
				end
			end
		end%if
	end%for
	
	% get type to export as
	exportTypeHandle = findobj('Tag','ExportAsType');
	exportType = get(exportTypeHandle, 'String');
	exportType = exportType{get(exportTypeHandle, 'Value')};
	for iC = clusterStart:length(MClust_Clusters)
		if ~isa(MClust_Clusters{iC}, exportType)
            MClust_Clusters{iC} = feval(exportType(2:end), GetName(MClust_Clusters{iC}), MClust_Clusters{iC});
		end
	end

	% close KKwik windows
	%CloseKK() - This is a pain.  Don't close. Feb 2008. ADR

	% go to Manual cutter if not open or update if open
	CutterFig = findobj('Type','figure','Tag', 'ClusterCutWindow');
   if isempty(CutterFig)
	   MClustCutter();
   else
	   MClustCutterClearClusterKeys(CutterFig);
	   MClustCutterRedrawClusterKeys(CutterFig);
   end
   
%%
case 'HoldCluster'
    
    ClusterUserData = get(cboHandle,'UserData');
    heldC = ClusterUserData{1};
    KlustaKwik_AllHold = ClusterUserData{2};
    
    % Check for holds before updating
    
    HoldNum = [];
    
    for iClust = 1:length(KlustaKwik_Clusters)
        HoldYN = get(findobj(figHandle, 'Tag', 'HoldCluster', 'UserData', {iClust, KlustaKwik_AllHold}), 'Value');
        if HoldYN % if we find a hold
            HoldNum = [HoldNum iClust];
        end
    end
    
	for iC = 1:length(KlustaKwik_Clusters)
        set(findobj(figHandle, 'Tag', 'HoldCluster', 'UserData', {iC, KlustaKwik_AllHold}), 'Value', 0);
	end
	if isempty(HoldNum)
		% we have turned off all holds
		% set correlations to off
		for iC = 1:length(KlustaKwik_Clusters)
			set(findobj(figHandle, 'Tag', 'Correlation', 'UserData', iC), 'String','----');
		end
	else
	 % turn self back on
        set(findobj(figHandle, 'Tag', 'HoldCluster', 'UserData', {heldC KlustaKwik_AllHold}), 'Value',1);
	% update correlations - added 3.5 ADR
	
        ClusterCorr = KKClust.WaveFormCorr(heldC,:);
		for iC = 1:length(KlustaKwik_Clusters)
			if iC ~= heldC
				set(findobj(figHandle, 'Tag', 'Correlation', 'UserData', iC), 'String',num2str(ClusterCorr(iC),'%.2f'));
			else
				set(findobj(figHandle, 'Tag', 'Correlation', 'UserData', iC), 'String','----');
			end
		end
	end
	
	KlustaKwikCallbacks('SelectCluster',cboHandle)

	% redraw if necessary
	CheckRedraw('Selection', figHandle);
    
case 'ChooseColor'
    col = get(cboHandle,'BackgroundColor');  % find the current color
    set(cboHandle,'backgroundcolor',get(0,'defaultuicontrolbackgroundcolor'))
    col = uisetcolor(col, 'Select a color'); % pass in the current color, if the user cancels, will not change color
    set(cboHandle, 'BackgroundColor', col);
	CheckRedraw('Selection', figHandle);

case 'SelectCluster'    
    % get ID
    iC = get(cboHandle, 'UserData');
    
    if ~isa(iC,'cell') % if the user data is a cell, this function has been called by HoldCluster
        set(cboHandle,'backgroundcolor',get(0,'defaultuicontrolbackgroundcolor'))
        selectedHandle = findobj(figHandle,'Tag','SelectCluster','BackgroundColor','c');
        if ~isempty(selectedHandle)
            set(selectedHandle,'backgroundcolor',get(0,'defaultuicontrolbackgroundcolor'))
        end
        set(cboHandle,'BackgroundColor','c')
        uicontrol(cboHandle);
    end
    
    CheckRedraw('Selection', figHandle);

    % Check for holds before updating window
    
    SelectedCluster = get(findobj(figHandle,'Tag','SelectCluster','BackgroundColor','c'),'UserData');
    
    HeldCluster = findobj(figHandle, 'Tag', 'HoldCluster', 'Value', 1);
    if isempty(HeldCluster)
        HoldNum = 0;
    else
        HoldNum = get(HeldCluster, 'UserData');
        HoldNum = HoldNum{1};
    end
       
    % change stats text
    objHandle = findobj(figHandle, 'Tag', 'StatisticsText');
    set(objHandle, 'String', KKClust.Stats{SelectedCluster});
    
    % plot average waveform
    AverageWaveformAxisHandle = findobj(figHandle, 'Tag', 'AverageWaveformAxis');
    axes(AverageWaveformAxisHandle); cla;
    
    if HoldNum
        mANDerr = KKClust.WaveForms{HoldNum};
        wfm = mANDerr{1};
        wferr = mANDerr{2};

		nWVSamples = size(wfm,2);
		
        for it = 1:4
            xrange = ((nWVSamples + 2) * (it-1)) + (1:nWVSamples); 
            hold on;
            plot(xrange, wfm(it,:));
            errorbar(xrange,wfm(it,:),wferr(it,:),'r'); 
        end
        hold on
    end
    
    mANDerr = KKClust.WaveForms{SelectedCluster};
    wfm = mANDerr{1};
    wferr = mANDerr{2};

	nWVSamples = size(wfm,2);
	
    for it = 1:4
        xrange = ((nWVSamples + 2) * (it-1)) + (1:nWVSamples); 
        hold on;
        plot(xrange, wfm(it,:));
        if MClust_ChannelValidity(it)
            errorbar(xrange,wfm(it,:),wferr(it,:)); 
        else
            errorbar(xrange,wfm(it,:),wferr(it,:),'k'); 
        end
    end

	set(gca,'Xlim',[0 4*(nWVSamples + 2)])
    title('Average Waveform');
    hold off
    set(gca, 'Tag', 'AverageWaveformAxis');
   
    % ISI histograms
    HistISIHandle = findobj(figHandle, 'Tag', 'ISIHistAxis');
    axes(HistISIHandle); cla;
    
    if HoldNum
        hist = KKClust.ISI{HoldNum};
        H = hist{1};
        binsUsed = hist{2};
        plot(binsUsed, H, 'r'); 
        hold on
    end
    
    hist = KKClust.ISI{SelectedCluster};
    H = hist{1};
    binsUsed = hist{2};
    plot(binsUsed, H); 
    xlabel('ISI, ms');
    set(gca, 'XScale', 'log');
    set(gca, 'YTick', max(H));
    set(gca, 'Tag', 'ISIHistAxis');
    hold off;   
    
  
case 'KeepCluster'
	CheckRedraw('Keep', figHandle);
	
% WINDOW PLOTTING FUNCTIONS
case {'xdim', 'ydim', 'Show0', 'ShowKeeps', 'CheckRedraw'}
	view2D = findobj(figHandle, 'Tag', 'View2D');
	if get(view2D, 'Value')
		Redraw2D(figHandle);
	end
	
	view3D = findobj(figHandle, 'Tag', 'View3D');
	if get(view3D, 'Value')
		Redraw3D(figHandle);
	end
	
	viewContour = findobj(figHandle, 'Tag', 'ViewContour');
	if get(viewContour, 'Value')
		RedrawContour(figHandle);
	end
	
case 'zdim'
	view3D = findobj(figHandle, 'Tag', 'View3D');
	if get(view3D, 'Value')
		Redraw3D(figHandle);
	end

case 'View2D'
	if get(cboHandle, 'Value')
		Redraw2D(figHandle);
	end
	
case 'View3D'
	RotateButton = findobj(figHandle, 'Tag', 'Rotate3D');
	if get(cboHandle, 'Value')
		set(RotateButton, 'Enable', 'on');
		Redraw3D(figHandle);
	else
		set(RotateButton, 'Enable', 'off');
	end		
	
case 'ViewContour'
	if get(cboHandle, 'Value')
		RedrawContour(figHandle);
	end
	
case 'Rotate3D'
    drawingFigHandle = findobj('Type', 'figure', 'Tag', 'KK3D');  % figure to draw in
    if ~isempty(drawingFigHandle)
        figure(drawingFigHandle);
        [az,el] = view;
        for iZ = 1:10:360
            figure(drawingFigHandle)
            view(az+iZ,el);
            drawnow;
        end
        set(cboHandle, 'value', 0);
    end
	
case 'Exit'
   close(figHandle);  
   
   drawingFigHandle = findobj('Type', 'figure', 'Tag', 'KK2D'); 
   if ~isempty(drawingFigHandle)
       close(drawingFigHandle);
   end    
   drawingFigHandle = findobj('Type', 'figure', 'Tag', 'KK3D');  
   if ~isempty(drawingFigHandle)
       close(drawingFigHandle);
   end    
   drawingFigHandle = findobj('Type', 'figure', 'Tag', 'KKContour');  
   if ~isempty(drawingFigHandle)
       close(drawingFigHandle);
   end

end


%------------------------------------
function CheckRedraw(key, figHandle)

view2D = get(findobj(figHandle, 'Tag', 'View2D'), 'Value');
view3D = get(findobj(figHandle, 'Tag', 'View3D'), 'Value');
viewContour = get(findobj(figHandle, 'Tag', 'ViewContour'), 'Value');
show0 = get(findobj(figHandle, 'Tag', 'Show0'), 'Value');
showKeeps = get(findobj(figHandle, 'Tag', 'ShowKeeps'), 'Value');
switch key
case 'Keep' % changed a keep
	if showKeeps && view2D; Redraw2D(figHandle); end
	if showKeeps && view3D; Redraw3D(figHandle); end
	if showKeeps && viewContour && ~show0; RedrawContour(figHandle); end	
case 'Selection' % only change selection
	if view2D; Redraw2D(figHandle); end
	if view3D; Redraw3D(figHandle); end
    if viewContour && ~show0; RedrawContour(figHandle); end
otherwise
    error('Reached otherwise in KlustaKwikCallbacks::CheckRedraw.');
end
	
%------------------------------------
function Redraw2D(figHandle)

% -- get variables
global KlustaKwik_Clusters  
global MClust_CurrentFeatures % used to keep track of which features are currently in memory
global MClust_CurrentFeatureNames % 
global MClust_CurrentFeatureData
global MClust_FeatureNames % names of features
global MClust_FeatureSources % <filenames, number pairs> for finding features in fd files

nClust = length(KlustaKwik_Clusters);

drawingFigHandle = findobj('Type', 'figure', 'Tag', 'KK2D');  % figure to draw in

global MClust_KK2D_Pos 

if isempty(drawingFigHandle)
	drawingFigHandle = figure('Name', '2D KlustaKwik viewer',...
        'NumberTitle', 'off', 'Tag', 'KK2D','Position', MClust_KK2D_Pos);
end

xdimHandle = findobj(figHandle, 'Tag', 'xdim'); xdim = get(xdimHandle, 'Value');           % x dimemsion to plot
ydimHandle = findobj(figHandle, 'Tag', 'ydim'); ydim = get(ydimHandle, 'Value');           % y dimension to plot
xlbls = get(xdimHandle, 'String');         % x labels
ylbls = get(ydimHandle, 'String');         % y lables

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

% figure already exists -- select it
figure(drawingFigHandle); cla; hold on;

% plot all points
if get(findobj(figHandle, 'Tag', 'Show0'), 'value')    
	h = plot(MClust_CurrentFeatureData(:,1), MClust_CurrentFeatureData(:,2), '.k','MarkerSize', 1);
end
% plot all keeps
if get(findobj(figHandle, 'Tag', 'ShowKeeps'), 'value')
	for iC = 1:length(KlustaKwik_Clusters)
		KeepYN = get(findobj(figHandle, 'Tag', 'KeepCluster', 'UserData', iC), 'Value');
		if KeepYN
			col = get(findobj(figHandle,'Tag','ChooseColor','UserData',iC),'Backgroundcolor');			
			[f,MCC] = FindInCluster(KlustaKwik_Clusters{iC});
			figure(drawingFigHandle);
			h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), '.', 'MarkerSize',1, 'Color', col);
		end
	end
end
% plot selection
selectedCluster = findobj(figHandle,'Tag','SelectCluster','Backgroundcolor', [0 1 1]);
if ~isempty(selectedCluster)
	iC = get(selectedCluster, 'UserData');
	col = get(findobj(figHandle,'Tag','ChooseColor','UserData',iC),'Backgroundcolor');			
	[f,MCC] = FindInCluster(KlustaKwik_Clusters{iC});
	figure(drawingFigHandle);
	h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), '.', 'MarkerSize',10, 'Color', col);
end
% plot hold
heldCluster = findobj(figHandle,'Tag','HoldCluster','Value',1);
if ~isempty(heldCluster)
	iC = get(heldCluster, 'UserData'); iC = iC{1};
	col = get(findobj(figHandle,'Tag','ChooseColor','UserData',iC),'Backgroundcolor');			
	[f,MCC] = FindInCluster(KlustaKwik_Clusters{iC});
	figure(drawingFigHandle);
	h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), 'x', 'MarkerSize',10, 'Color', col);
end
figure(drawingFigHandle);
set(gca, 'XLim', [min(MClust_CurrentFeatureData(:,1)) max(MClust_CurrentFeatureData(:, 1))+0.0001]);
set(gca, 'YLim', [min(MClust_CurrentFeatureData(:,2)) max(MClust_CurrentFeatureData(:, 2))+0.0001]);
xlabel(MClust_CurrentFeatureNames{1}, 'interpreter', 'none');
ylabel(MClust_CurrentFeatureNames{2}, 'interpreter', 'none');
zoom on

%-------------------------------------------------------
function Redraw3D(figHandle)

% -- get variables
global KlustaKwik_Clusters 
global MClust_CurrentFeatures % used to keep track of which features are currently in memory
global MClust_CurrentFeatureNames % 
global MClust_CurrentFeatureData
global MClust_FeatureNames % names of features
global MClust_FeatureSources % <filenames, number pairs> for finding features in fd files

nClust = length(KlustaKwik_Clusters);

global MClust_KK3D_Pos

drawingFigHandle = findobj('Type', 'figure', 'Tag', 'KK3D');  % figure to draw in
if isempty(drawingFigHandle)
	drawingFigHandle = figure('Name', '3D KlustaKwik viewer',...
        'NumberTitle', 'off', 'Tag', 'KK3D','Position', MClust_KK3D_Pos);
end

xdimHandle = findobj(figHandle, 'Tag', 'xdim'); xdim = get(xdimHandle, 'Value');           % x dimemsion to plot
ydimHandle = findobj(figHandle, 'Tag', 'ydim'); ydim = get(ydimHandle, 'Value');           % y dimension to plot
zdimHandle = findobj(figHandle, 'Tag', 'zdim'); zdim = get(zdimHandle, 'Value');           % y dimension to plot
xlbls = get(xdimHandle, 'String');         % x labels
ylbls = get(ydimHandle, 'String');         % y labels
zlbls = get(zdimHandle, 'String');         % z labels

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
% get ydim
if (MClust_CurrentFeatures(3) ~= zdim)
    temp = load(MClust_FeatureSources{ydim,1}, '-mat', 'FeatureData');
    MClust_CurrentFeatureData(:,3) = temp.FeatureData(:,MClust_FeatureSources{zdim,2});
    MClust_CurrentFeatures(3) = zdim;
    MClust_CurrentFeatureNames{3} = MClust_FeatureNames{zdim};
end

% figure already exists -- select it
figure(drawingFigHandle); cla; hold on;

% plot all points
if get(findobj(figHandle, 'Tag', 'Show0'), 'value')
	h = plot3(MClust_CurrentFeatureData(:,1), MClust_CurrentFeatureData(:,2), MClust_CurrentFeatureData(:,3), '.k','MarkerSize', 1);
end
% plot all keeps
if get(findobj(figHandle, 'Tag', 'ShowKeeps'), 'value')
	for iC = 1:length(KlustaKwik_Clusters)
		KeepYN = get(findobj(figHandle, 'Tag', 'KeepCluster', 'UserData', iC), 'Value');
		if KeepYN
			col = get(findobj(figHandle,'Tag','ChooseColor','UserData',iC),'Backgroundcolor');			
			[f,MCC] = FindInCluster(KlustaKwik_Clusters{iC});
			figure(drawingFigHandle);
			h = plot3(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), MClust_CurrentFeatureData(f,3), '.', 'MarkerSize',1, 'Color', col);
		end
	end
end
% plot selection
selectedCluster = findobj(figHandle,'Tag','SelectCluster','Backgroundcolor', [0 1 1]);
if ~isempty(selectedCluster)
	iC = get(selectedCluster, 'UserData');
	col = get(findobj(figHandle,'Tag','ChooseColor','UserData',iC),'Backgroundcolor');			
	[f,MCC] = FindInCluster(KlustaKwik_Clusters{iC});
	figure(drawingFigHandle);
	h = plot3(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), MClust_CurrentFeatureData(f,3), '.', 'MarkerSize',10, 'Color', col);
end
% plot hold
heldCluster = findobj(figHandle,'Tag','HoldCluster','Value',1);
if ~isempty(heldCluster)
	iC = get(heldCluster, 'UserData'); iC = iC{1};
	col = get(findobj(figHandle,'Tag','ChooseColor','UserData',iC),'Backgroundcolor');			
	[f,MCC] = FindInCluster(KlustaKwik_Clusters{iC});
	figure(drawingFigHandle);
	h = plot3(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), MClust_CurrentFeatureData(f,3), 'x', 'MarkerSize',10, 'Color', col);
end
figure(drawingFigHandle);
set(gca, 'XLim', [min(MClust_CurrentFeatureData(:,1)) max(MClust_CurrentFeatureData(:, 1))+0.0001]);
set(gca, 'YLim', [min(MClust_CurrentFeatureData(:,2)) max(MClust_CurrentFeatureData(:, 2))+0.0001]);
set(gca, 'ZLim', [min(MClust_CurrentFeatureData(:,3)) max(MClust_CurrentFeatureData(:, 3))+0.0001]);
xlabel(MClust_FeatureNames{xdim});
ylabel(MClust_FeatureNames{ydim});
zlabel(MClust_FeatureNames{zdim});
view(3)
rotate3D on

%--------------------------------------------------
function RedrawContour(figHandle)

% -- get variables
global KlustaKwik_Clusters MClust_FeatureData 
global MClust_CurrentFeatures % used to keep track of which features are currently in memory
global MClust_CurrentFeatureNames % 
global MClust_CurrentFeatureData
global MClust_FeatureNames % names of features
global MClust_FeatureSources % <filenames, number pairs> for finding features in fd files


nClust = length(KlustaKwik_Clusters);

global MClust_KKContour_Pos

drawingFigHandle = findobj('Type', 'figure', 'Tag', 'KKContour');  % figure to draw in
if isempty(drawingFigHandle)
	drawingFigHandle = figure('Name', 'Contour KlustaKwik viewer',...
        'NumberTitle', 'off', 'Tag', 'KKContour','Position', MClust_KKContour_Pos);
end

xdimHandle = findobj(figHandle, 'Tag', 'xdim'); xdim = get(xdimHandle, 'Value');           % x dimemsion to plot
ydimHandle = findobj(figHandle, 'Tag', 'ydim'); ydim = get(ydimHandle, 'Value');           % y dimension to plot
xlbls = get(xdimHandle, 'String');         % x labels
ylbls = get(ydimHandle, 'String');         % y lables

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

% figure already exists -- select it
figure(drawingFigHandle); cla;

XData = []; 
YData = [];
% plot all points 
if get(findobj(figHandle, 'Tag', 'Show0'), 'value')
	XData = MClust_CurrentFeatureData(:,1)';
	YData = MClust_CurrentFeatureData(:,2)';
else
	f = [];
	% for all clusters include y/n
	% keeps
	if get(findobj(figHandle, 'Tag', 'ShowKeeps'), 'value');
		for iC = 1:nClust
			if get(findobj(figHandle, 'Tag', 'KeepCluster', 'UserData', iC), 'value')
				f = cat(1, f, FindInCluster(KlustaKwik_Clusters{iC}));
			end
		end
	end
	selectedCluster = findobj(figHandle,'Tag','SelectCluster','Backgroundcolor', [0 1 1]);
	if ~isempty(selectedCluster)
		iC = get(selectedCluster, 'UserData');
		f = cat(1, f, FindInCluster(KlustaKwik_Clusters{iC}));
	end
	heldCluster = findobj(figHandle,'Tag','HoldCluster','Value',1);
	if ~isempty(heldCluster)
		iC = get(heldCluster, 'UserData'); iC = iC{1};
		f = cat(1, f, FindInCluster(KlustaKwik_Clusters{iC}));
	end
    if ~isempty(f)
        XData = MClust_CurrentFeatureData(f,1)';
        YData = MClust_CurrentFeatureData(f,2)';
    end
end
% make the plot
nX = 100; nY = 100;
H = [];

xmin = min(MClust_CurrentFeatureData(:,1))+0.0001; 
xmax = max(MClust_CurrentFeatureData(:,1))+0.0001; 
ymin = min(MClust_CurrentFeatureData(:,2))+0.0001; 
ymax = max(MClust_CurrentFeatureData(:,2))+0.0001; 

if ~isempty(XData)
    H = ndhist([XData; YData], [nX; nY], [xmin; ymin], [xmax; ymax])';
    H0 = Hsmooth(H);
    H0 = fixedLog(H,-1);
    if any(H0(:)>-1)
        H = H0;
    end
    [nX, nY] = size(H);
    X = linspace(xmin, xmax,nX);
    Y = linspace(ymin, ymax,nY);
end

figure(drawingFigHandle);
if isempty(XData) || isempty(H) || all(H(:)==H(1))    
    cla;
    set(gca, 'XTick', [], 'YTick', []);
else
    Cont = contourf(X, Y , H, 15);
    set(gca, 'xlim', [xmin xmax], 'ylim', [ymin ymax], 'XTick', [], 'YTick', []);
end
xlabel(MClust_FeatureNames{xdim});
ylabel(MClust_FeatureNames{ydim});
zoom on;

%==========================================================
function fL = fixedLog(H, infimum)

%  returns the log(H) with H a N-dim array 
%  where the NaNs of log(x<=0) of H are replaced by infimum
%

fL = H;
fL(fL <0) = 0;
fL = log(fL);
fL(fL <= infimum) = infimum;
return;


%==========================================
function [S, bb] = Hsmooth(H)
%
%  [S, bb] = Hsmooth(H)
%
% smoothes a 2D histogram H (= 2D array)
% with a 6-th order double low pass firls bb (linear phase least-square FIR filter)
%

% create filter
if ~isempty(which('firls'))
	b = firls(6, [0 .5 .5 1], [1 .5 .5 0]);  
	bb = kron(b',b);    % 2D filter = tensor product of 1D filters
	
	S = filter2(bb,H);  % first pass (introduces a linear phase shift)
	S = filter2(bb',S);  % second pass (compensates phase shift)
else
	S = H;
end

%===========================================
function CloseKK()

drawingFigHandle2D = findobj('Type', 'figure', 'Tag', 'KK2D');
if ~isempty(drawingFigHandle2D)
	close(drawingFigHandle2D);
end

drawingFigHandle3D = findobj('Type', 'figure', 'Tag', 'KK3D');
if ~isempty(drawingFigHandle3D)
	close(drawingFigHandle3D);
end

drawingFigHandleContour = findobj('Tag', 'KKContour');
if ~isempty(drawingFigHandleContour)
	close(drawingFigHandleContour);
end

KKDW = findobj('Tag', 'KKDecisionWindow');
if ~isempty(KKDW)
	close(KKDW);
end
