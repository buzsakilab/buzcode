function MClustCallbacks(varargin)

% MClustCallbacks
% 
% Contains all callbacks for the BBClust main window (see BBClust)
%
% This program stores a number of key parameters as global variables.
% All global variables start with the tag "MClust_".  
% 
% See also
%    MClust
%
% ADR 1998, 2007
% ADR 2007
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.
% ADR Jan/2007

% global variables
global MClust_Directory     % directory where the MClust main function and features reside
global MClust_TTfn         % file name for tt file
global MClust_FDfn         % file name for fd file
global MClust_TTdn         % directory name for tt file
global MClust_FDdn         % directory name for fd file
global MClust_TText        % extension for tt file
global MClust_FDext        % extension for fd file
global MClust_CurrentFeatures % used to keep track of which features are currently in memory
global MClust_FeatureData  % features calculated from tt data
global MClust_FeatureSources
global MClust_FeatureNames % names of features
global MClust_FeaturesToUse % names of features to use (for feature calculation)
global MClust_FeatureTimestamps 
global MClust_ChannelValidity % 4 x 1 array of channel on (1) or off (0) flags
global MClust_Clusters     % cell array of cluster objects 
global MClust_FilesWrittenYN
global MClust_NeuralLoadingFunction

global KlustaKwik_Clusters

%% find ssw key
if ~isempty(varargin)
  ssw = varargin(1);
  ssw = ssw{1};
  MClustFigureHandle = findobj('Tag', 'MClustMainWindow');
  cboHandle = [];
else
  cboHandle = gcbo;                          % current uicontrol handle
  MClustFigureHandle = ParentFigureHandle(cboHandle); % current figure handle

  ssw = get(cboHandle, 'Tag');
end

switch ssw
%% Select Loading Engine
	case 'SelectLoadingEngine'
		%requires MClust_NeuralLoadingFunction
		value = get(cboHandle, 'Value');
		LoadingEngines = get(cboHandle, 'String');
		MClust_NeuralLoadingFunction = LoadingEngines{value};
	
%% Load Features
	case 'LoadFeaturesButton'

		% initialize
		set(findobj(MClustFigureHandle, 'Tag', 'LoadFeaturesButton'), 'Value', 0);
		
		% get features to use
		featuresToUse = GetFeaturesFromListBox(MClustFigureHandle);
		
		%find basename
		[fn dn] = uigetfile('*T*.dat;*.ntt;*.nse;*.nst;*.tt', ...
			'Select the spike data file from the desired tetrode');
		if fn == 0    % User hit cancel
			return
		end
		if ~isempty(dn)
			dn = dn(1:end-1);  % remove filesep
		end
		MClust_TTdn = dn;
		[p MClust_TTfn MClust_TText] = fileparts(fn);       	

		MClust_FDfn = [MClust_TTfn MClust_FDext];
        pushdir(MClust_TTdn);
        if exist('FD','dir')
            MClust_FDdn = [dn filesep 'FD'];
        else
            MClust_FDdn = dn;
        end
        popdir;
  
		% Recalc
        TTFileNameHandle = findobj(MClustFigureHandle, 'Tag', 'TTFileName');
		set(TTFileNameHandle,'BackGroundColor','r');
		if CalculateFeatures(MClust_TTfn, featuresToUse) % calculateOnlyIfNecessary
			set(TTFileNameHandle,'BackGroundColor','c');
			set(findobj(MClustFigureHandle, 'Tag', 'LoadFeaturesButton'), 'Value', 1);
        end
        
        % put filename in textfile
        set(TTFileNameHandle, 'string', MClust_TTfn);
		

%% Transfer between features to use and features to ignore
	case {'FeaturesUseListbox', 'FeaturesIgnoreListbox'}
		TransferBetweenListboxes;
		FeaturesUseListbox = findobj(MClustFigureHandle, 'Tag', 'FeaturesUseListbox');		
		ChooseFeaturesButton = findobj(MClustFigureHandle, 'Tag', 'ChooseFeaturesButton');       
        
		featuresToUse = GetFeaturesFromListBox(MClustFigureHandle);

		if ~isempty(MClust_TTfn)
			% we have a TT in progress -- need to recalc
			set(findobj(MClustFigureHandle, 'Tag', 'TTFileName'),'BackGroundColor','r');
			if CalculateFeatures(MClust_TTfn, featuresToUse) % calculateOnlyIfNecessary
				set(findobj(MClustFigureHandle, 'Tag', 'TTFileName'),'BackGroundColor','c');
				set(findobj(MClustFigureHandle, 'Tag', 'LoadFeaturesButton'), 'Value', 1);
			end
		end

			
%% Channel Validity
	case 'TTValidity1'
		MClust_ChannelValidity(1) = get(cboHandle, 'Value');

	case 'TTValidity2'
		MClust_ChannelValidity(2) = get(cboHandle, 'Value');

	case 'TTValidity3'
		MClust_ChannelValidity(3) = get(cboHandle, 'Value');

	case 'TTValidity4'
		MClust_ChannelValidity(4) = get(cboHandle, 'Value');

%% KlustaKwik
case 'KlustaKwikSelection'
	
	if isempty(MClust_FeatureNames)
		errordlg('No features calculated.', 'MClust error', 'modal');
		return
	end
	
	% check for klustakwik window
	KKwin = findobj('Tag', 'KKDecisionWindow');
	if ~isempty(KKwin)
		warndlg('KlustaKwik window still open.', 'KlustaKwik/MClust', 'modal');
		return
	end
	% else
    
    pushdir(MClust_FDdn);
    expExtension = [MClust_TTfn '.clu.1'];
    CLUfiles = FindFiles([MClust_TTfn '.clu.1']); % Can we figure out the correct file?
    if length(CLUfiles) ~= 1     % NO
        [fn, fdir] = uigetfile(expExtension);
        % Get the file and remove any . extensions (there are 
        % more than one so you need to run this a couple times).
    else                         % YES
        [fdir fn fext] = fileparts(CLUfiles{1});
    end
    popdir;
    
    if isempty(fn)
        return
    end
    
    [p rootname e ] = fileparts(fn);
    [p rootname e ] = fileparts(rootname);
    [p rootname e ] = fileparts(rootname);
    
    
    if rootname       
		% Load in the clusters from KlustaKwik
		file_no = 1;
		clu_file = [fullfile(fdir,rootname) '.clu.' num2str(file_no)];
		KlustaKwik_Clusters = KlustaKwikImport(clu_file);
    end
       
    if ~isempty(KlustaKwik_Clusters)
        KlustaKwikDecisionWindow; % create window
    end

%% Generalized Cutter
    case 'CutPreClusters'
		if isempty(MClust_FeatureNames)
			errordlg('No features calculated.', 'MClust error', 'modal');
		else
			MClustCutter;
		end
		
%% Save defaults
	case 'SaveDefaults'
		[fn,dn] = uiputfile(fullfile(MClust_Directory,'defaults.mclust'));
		if fn, SaveDefaults(MClustFigureHandle, dn,fn); end

%% Load defaults
	case 'LoadDefaults'
		[fn,dn] = uigetfile('*.mclust');
		if fn, LoadDefaults(MClustFigureHandle, fullfile(dn,fn)); end
		
%% Clear Clusters
	case 'ClearClusters'
		
		ynClose = questdlg('Clearing clusters.  No undo available. Are you sure?', 'ClearQuestion', 'Yes', 'Cancel', 'Cancel');
		if strcmp(ynClose,'Yes')
			ClusterCutWindowHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');
			if ~isempty(ClusterCutWindowHandle)
				MClustCutterClearClusterKeys(ClusterCutWindowHandle)
			end
			
			MClust_Clusters = {};			

			if ~isempty(ClusterCutWindowHandle)
				MClustCutterRedrawClusterKeys(ClusterCutWindowHandle)
				
				redrawAxesHandle = findobj(ClusterCutWindowHandle, 'Tag', 'RedrawAxes');
				redrawAxesFlag = get(redrawAxesHandle, 'Value');
				if redrawAxesFlag
					MClustCutterCallbacks('RedrawAxes');
				end
			end
		end
		
%% Save Clusters
	case 'SaveClusters'
		OK = SaveClusters();
		if OK
			msgbox('Clusters saved successfully.', 'MClust msg');
		else
			errordlg('Clusters not saved.', 'MClust error', 'modal');
		end
		
%% Load Clusters
	case 'LoadClusters'
		if isempty(MClust_TTfn)
			errordlg('Select a TT and calculate features first before loading clusters.',...
				'MClust error', 'modal');
		else
			OK = LoadClusters();
			if OK
				msgbox('Clusters appended to clusters list.', 'MClust msg');
			else
				errordlg('Clusters not loaded.', 'MClust error', 'modal');
			end
        end

%% Load Clusters
	case 'ApplyConvexHulls'
        MClust_ApplyConvexHulls;
        ClusterCutWindowHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');
        if ~isempty(ClusterCutWindowHandle)
            MClustCutterClearClusterKeys(ClusterCutWindowHandle)
            MClustCutterRedrawClusterKeys(ClusterCutWindowHandle)

            redrawAxesHandle = findobj(ClusterCutWindowHandle, 'Tag', 'RedrawAxes');
            redrawAxesFlag = get(redrawAxesHandle, 'Value');
            if redrawAxesFlag
                MClustCutterCallbacks('RedrawAxes');
            end
        end

%% Clear Workspace
	case 'ClearWorkspaceOnly'
		if ~strcmp(MClust_FilesWrittenYN,'yes')
			ynWrite = questdlg('Are you sure?  No undo is possible', 'ExitQuestion', 'Yes', 'Cancel','Cancel');
			if streq(ynWrite, 'Yes'),
				ClearWorkspace();
			end
		else
			ClearWorkspace();
		end

%% Exit

	case 'ExitOnlyButton'
		% requires MClust_FilesWrittenYN

		if ~strcmp(MClust_FilesWrittenYN,'yes')
			ynWrite = questdlg('Are you sure?', 'ExitQuestion', 'Yes','Cancel', 'Cancel');
			if streq(ynWrite, 'Cancel'), return; end
		end
		ClearWorkspace;
		clear global MClust_* MCLUST_*
		close(MClustFigureHandle);

%% WriteFiles
	case 'WriteFiles'
		if isempty(MClust_Clusters)
			errordlg('No clusters exist.', 'MClust error', 'modal');
			return
		end
		SaveClusters([fullfile(MClust_TTdn, MClust_TTfn), '.clusters']);
		OK = WriteTFiles(MClust_Clusters);
		if OK 
			msgbox('Files written.', 'MClust msg');
			MClust_FilesWrittenYN = 'yes';
		else
			msgbox('Files NOT written.', 'MClust msg');
		end

%%
	otherwise
		warndlg('Sorry, feature not yet implemented.');
end % switch

end

%% SaveClusters
function OK = SaveClusters(fn)
global MClust_Clusters MClust_Colors 
global MClust_TTdn MClust_TTfn

if nargin==0
	[fn,dn] = uiputfile(fullfile(MClust_TTdn, [MClust_TTfn '.clusters']));
	fn = fullfile(dn, fn);
end
if fn
	save(fn, 'MClust_Clusters','MClust_Colors','-mat');
	OK = true;
else
	OK = false;
end

end

%% LoadClusters
function OK = LoadClusters(fn)
global MClust_Clusters MClust_Colors 
global MClust_TTdn MClust_TTfn

if nargin==0
	[fn,dn] = uigetfile(fullfile(MClust_TTdn, [MClust_TTfn '.clusters']));
	fn = fullfile(dn, fn);
end
if fn
	temp = load(fn,'-mat');
	nClusters = length(MClust_Clusters);
	nToAdd = length(temp.MClust_Clusters);
	MClust_Clusters = cat(2, MClust_Clusters, temp.MClust_Clusters);
	MClust_Colors((nClusters+2):(nClusters+nToAdd+1),:) = temp.MClust_Colors(2:(nToAdd+1),:);

    ClusterCutWindow = findobj('Type','figure','Tag', 'ClusterCutWindow');
    if ~isempty(ClusterCutWindow)
        close(ClusterCutWindow);
    end
  
	CHDrawingAxisWindow = findobj('Type','figure','Tag', 'CHDrawingAxisWindow');
	if ~isempty(CHDrawingAxisWindow)
		close(CHDrawingAxisWindow);
	end
	
	OK = true;
else
	OK = false;
end

end

%% SaveDefaults
function SaveDefaults(MClustFigureHandle, dn,fn)
	
   global MClust_Colors MClust_ChannelValidity MClust_Directory MClust_ClusterSeparationFeatures
   global MClust_ClusterCutWindow_Pos MClust_CHDrawingAxisWindow_Pos MClust_KKDecisionWindow_Pos
   global MClust_KK2D_Pos MClust_KK3D_Pos MClust_KKContour_Pos MClust_ClusterCutWindow_Marker
   global MClust_ClusterCutWindow_MarkerSize MClust_NeuralLoadingFunction
   
   featuresToUse = GetFeaturesFromListBox(MClustFigureHandle);

   save(fullfile(dn,fn), ...
          'MClust_Colors', ...
          'featuresToUse', 'MClust_ChannelValidity',  ...
          'MClust_NeuralLoadingFunction','MClust_ClusterSeparationFeatures', ...
          'MClust_ClusterCutWindow_Pos','MClust_CHDrawingAxisWindow_Pos',...
          'MClust_KKDecisionWindow_Pos','MClust_KK2D_Pos','MClust_KK3D_Pos','MClust_KKContour_Pos','MClust_ClusterCutWindow_Marker', ...
          'MClust_ClusterCutWindow_MarkerSize','-mat');
	  msgbox('Defaults saved successfully.', 'MClust msg');
end

%% LoadDefaults
function LoadDefaults(MClustFigureHandle, fn)
   global MClust_Colors MClust_ChannelValidity MClust_Directory MClust_ClusterSeparationFeatures
   global MClust_ClusterCutWindow_Pos MClust_CHDrawingAxisWindow_Pos MClust_KKDecisionWindow_Pos
   global MClust_KK2D_Pos MClust_KK3D_Pos MClust_KKContour_Pos MClust_ClusterCutWindow_Marker
   global MClust_ClusterCutWindow_MarkerSize MClust_NeuralLoadingFunction
   load(fn, '-mat');

   % fix features to use boxen
   uifeaturesIgnoreHandle = findobj(MClustFigureHandle, 'Tag', 'FeaturesIgnoreListbox');
   uifeaturesUseHandle = findobj(MClustFigureHandle, 'Tag', 'FeaturesUseListbox');
   uiChooseFeaturesButton = findobj(MClustFigureHandle, 'Tag', 'ChooseFeatures');
   allFeatures = [get(uifeaturesIgnoreHandle, 'String'); get(uifeaturesUseHandle, 'String')];
   featureIgnoreString = {};
   featureUseString = {};
   for iF = 1:length(allFeatures)
	   if any(strcmp(allFeatures{iF}, featuresToUse))
		   featureUseString = cat(1, featureUseString, allFeatures(iF));
	   else
		   featureIgnoreString = cat(1, featureIgnoreString, allFeatures(iF));
	   end
   end
   set(uifeaturesIgnoreHandle, 'String', featureIgnoreString);
   set(uifeaturesUseHandle, 'String', featureUseString);
   if ~isempty(featureUseString)
	   set(uiChooseFeaturesButton, 'Value', 1);
   end
   
   % neural loading function
   uiLoadingEngine = findobj(MClustFigureHandle, 'Tag', 'SelectLoadingEngine');
   value = strmatch(MClust_NeuralLoadingFunction, get(uiLoadingEngine, 'String'));
   set(uiLoadingEngine, 'value', value);

   % fix channel validity
   set(findobj(MClustFigureHandle, 'Tag', 'TTValidity1'), 'value', MClust_ChannelValidity(1));
   set(findobj(MClustFigureHandle, 'Tag', 'TTValidity2'), 'value', MClust_ChannelValidity(2));
   set(findobj(MClustFigureHandle, 'Tag', 'TTValidity3'), 'value', MClust_ChannelValidity(3));
   set(findobj(MClustFigureHandle, 'Tag', 'TTValidity4'), 'value', MClust_ChannelValidity(4));
   
   msgbox('Defaults loaded successfully.', 'MClust msg');

 end

%% function Clear Workspace
function ClearWorkspace()
    global MClust_TTdn 
    
	MClustFigureHandle = findobj('Tag', 'MClustMainWindow');
    
    % close all windows
    KKWindow = findobj('Type', 'figure', 'Tag', 'KKDecisionWindow');
    if ~isempty(KKWindow)
        KlustaKwikCallbacks('Exit');
    end
    ClusterCutWindow = findobj('Type','figure','Tag', 'ClusterCutWindow');
    if ~isempty(ClusterCutWindow)
        MClustCutterCallbacks('Exit');
    end

    if ~isempty(MClust_TTdn)
        popdir all;
        pushdir(MClust_TTdn);
    end
    
    % reset globals
	MClustResetGlobals('workspace', MClustFigureHandle);
	
	set(findobj(MClustFigureHandle, 'Tag', 'LoadFeaturesButton'), 'Value', 0,'Enable','on');
	set(findobj(MClustFigureHandle, 'Tag', 'FeaturesIgnoreListbox'), 'Enable','on');
	set(findobj(MClustFigureHandle, 'Tag', 'FeaturesUseListbox'), 'Enable','on');
    
    % update Channel Validity Checkboxes with used channels in featurefile
    set(findobj(MClustFigureHandle, 'Tag', 'TTValidity1'), 'Value', 1, 'Enable', 'on');
    set(findobj(MClustFigureHandle, 'Tag', 'TTValidity2'), 'Value', 1, 'Enable', 'on');
    set(findobj(MClustFigureHandle, 'Tag', 'TTValidity3'), 'Value', 1, 'Enable', 'on');
    set(findobj(MClustFigureHandle, 'Tag', 'TTValidity4'), 'Value', 1, 'Enable', 'on');
    
    set(findobj(MClustFigureHandle, 'Tag', 'TTFileName'), 'String', [],'BackGroundColor',[0.7 0.7 0.7]);
    
end

%% function GetFeaturesFromListbox
function featuresToUse = GetFeaturesFromListBox(MClustFigureHandle)
      featureToUseHandle = findobj(MClustFigureHandle, 'Tag', 'FeaturesUseListbox');
      featuresToUse = get(featureToUseHandle, 'String');
end
