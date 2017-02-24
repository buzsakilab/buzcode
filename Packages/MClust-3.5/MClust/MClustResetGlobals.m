function MClustResetGlobals(flag, MClustFigureHandle)

% MClustResetGlobals('workspace', MClustFigureHandle)
% MClustResetGlobals('initialize');
%
% initialize sets ALL MClust globals to the initial values
% workspace clears only globals involved in the workspace
%
% ADR MClust 3.5 Feb 2008

% full list of globals used in MClust
% -- constant globals
global MCLUST_DEBUG         % if true then prints verbose debugging comments
global MCLUST_VERSION       % version flag

% -- globals initialized once
global MClust_Directory     % directory where the MClust main function and features reside
global MClust_AvailableClusterTypes

% -- globals used to communicate between functions
global MClust_FeatureData  % feature data
global MClust_TTfn         % file name for tt file
global MClust_FDfn         % file name for fd file
global MClust_TTdn         % directory name for tt file
global MClust_FDdn         % directory name for fd file
global MClust_TText        % extension for tt file
global MClust_FDext        % extension for fd file

global MClust_Clusters     % the cluster objects
global MClust_Colors       % colors for each cluster 
global MClust_Hide         % Hide or Show (1/0) for each cluster
global MClust_UnaccountedForOnly % is the 0 cluster showing everything?
global MClust_FilesWrittenYN % have the files been written

global MClust_Undo MClust_Redo

global MClust_ClusterSeparationFeatures % which features is CQ being run on?
global MClust_max_records_to_load

% window positions and limits
global MClust_ClusterCutWindow_Pos
global MClust_CHDrawingAxisWindow_Pos
global MClust_KKDecisionWindow_Pos
global MClust_KK2D_Pos 
global MClust_KK3D_Pos
global MClust_KKContour_Pos
global MClust_ClusterCutWindow_Marker
global MClust_ClusterCutWindow_MarkerSize
global MClust_AverageWaveform_ylim;

global MClust_TTData       % data from tt file
global MClust_CurrentFeatures % used to keep track of which features are currently in memory
global MClust_FeatureNames % names of features
global MClust_FeatureSources % <filenames, number pairs> for finding features in fd files
global MClust_FeaturesToUse % features to use
global MClust_FeatureTimestamps % timestamps
global MClust_ChannelValidity % 4 x 1 array of channel on (1) or off (0) flags
global MClust_NeuralLoadingFunction % Loading Engine 

global KlustaKwik_Clusters % clustes being used by KKwik
global KKC_IDSet
global KKClust

switch flag
	case 'initialize'
		% constants
		MCLUST_DEBUG = true;
		MCLUST_VERSION = '3.5A.05, 28 July 2008';

		MClust_Directory = fileparts(which('MClust.m'));
		
		% Find available Cluster Types
		ClusterTypeDirs = dir(fullfile(MClust_Directory, 'ClusterTypes\@*'));
		MClust_AvailableClusterTypes = cell(length(ClusterTypeDirs),1);
		for iC = 1:length(ClusterTypeDirs)
			MClust_AvailableClusterTypes{iC} = ClusterTypeDirs(iC).name(2:end);
		end

		% defaults
		MClust_TText = '.ntt';
		MClust_FDext = '.fd';
        
        MClust_FeatureData = []; % feature data
        MClust_TTfn = '';       % file name for tt file
        MClust_FDfn = '';       % file name for fd file
        MClust_TTdn = '';      % directory name for tt file
        MClust_FDdn = '';        % directory name for fd file

		MClust_Clusters = {};

		MClust_Colors = colorcube; close;

		MClust_Hide = zeros(99,1);
		MClust_UnaccountedForOnly = 0;

		MClust_FilesWrittenYN = 'no';

		MClust_ClusterSeparationFeatures = {'energy','wavePC1'};

		MClust_ClusterCutWindow_Pos= [10 60 450 650];
		MClust_CHDrawingAxisWindow_Pos= [500 200 650 650];
		MClust_KKDecisionWindow_Pos = [10 60 600 780];
		MClust_KK2D_Pos = [629 379 515 409];
		MClust_KK3D_Pos = [633 300 510 363];
		MClust_KKContour_Pos = [702 44 391 249];

		MClust_ClusterCutWindow_Marker = 1;
		MClust_ClusterCutWindow_MarkerSize = 1;
		
		MClust_Undo = {}; MClust_Redo = {};

		MClust_AverageWaveform_ylim = [-2100 2100];
        
        MClust_max_records_to_load = 10000; % ADR 26 Feb 2008

        MClust_TTData = [];      % data from tt file
        MClust_CurrentFeatures = [-1,-1,-1];% used to keep track of which features are currently in memory
        MClust_FeatureNames = {}; % names of features
        MClust_FeatureSources = {}; % <filenames, number pairs> for finding features in fd files
        MClust_FeatureTimestamps = []; % timestamps
        MClust_ChannelValidity = [1 1 1 1];% 4 x 1 array of channel on (1) or off (0) flags
        MClust_NeuralLoadingFunction = []; % Loading Engine

        KlustaKwik_Clusters = {};% clustes being used by KKwik
        KKC_IDSet = {}; KKClust = {};

	case 'workspace'
		% modified to save variables needed after clearing - 3.5 ncst

		clear global KlustaKwik_* MClust_Fe* MClust_TTD* MClust_TTfn
		clear global MClust_CurrentFeatureData
		clear global KKC_IDSet KKClust 
        clear global GC featureindex featuresToUse

        MClust_FilesWrittenYN = 'no';

		MClust_Clusters = {};

		MClust_Hide = zeros(99,1);
		MClust_UnaccountedForOnly = 0;
		
		MClust_Undo = {}; MClust_Redo = {};

		FeaturesUseListbox = findobj(MClustFigureHandle, 'Tag', 'FeaturesUseListbox');
		FeaturesToUse = get(FeaturesUseListbox, 'String')';
		MClust_FeaturesToUse = FeaturesToUse;
        MClust_max_records_to_load = 10000; % ADR Feb 2008
        
        MClust_TTData = [];      % data from tt file
        MClust_CurrentFeatures = [-1,-1,-1];% used to keep track of which features are currently in memory
        MClust_FeatureNames = {}; % names of features
        MClust_FeatureSources = {}; % <filenames, number pairs> for finding features in fd files
        MClust_FeatureTimestamps = []; % timestamps
        MClust_ChannelValidity = [1 1 1 1];% 4 x 1 array of channel on (1) or off (0) flags
      
end
