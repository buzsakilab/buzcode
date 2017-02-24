function ok = CalculateFeatures(basename, featureList)

% CalculateFeatures(basename, featureList)
%
% Goes through featureList and makes sure features are calculated for
% electrode file basename.
%
% Recalculates them if necessary
%
% ADR 2007
% Version M3.5.
% Released with MClust 3.5.

% PARAMETERS
global MClust_FDdn MClust_FDext MClust_ChannelValidity MCLUST_DEBUG
global MClust_max_records_to_load MClust_TTdn MClust_TTfn MClust_TText
global MClust_FeatureTimestamps

global MClust_FeatureNames % names of features
global MClust_FeatureSources % <filenames, number pairs> for finding features in fd files

if ~isempty(MClust_max_records_to_load)
	record_block_size = MClust_max_records_to_load;
else
	record_block_size = 30000;  % maximum number of spikes to load into memory at one time
end

template_matching = 0; % used to remove noise spikes which are not "spike-like,", template-matching is not a currently supported function

NormalizeFDYN = 'no'; % 'yes' if you want to normalize the feature data files to mean = 0, std = 1

if isempty(featureList)
	error('MClust:InternalError', 'Empty feature list sent to CalculateFeatures.');
end

% GO

featureFiles = cell(size(featureList));
for iF = 1:length(featureList)
	featureFiles{iF} = fullfile(MClust_FDdn, [basename '_' featureList{iF}]);
end

if MCLUST_DEBUG	
	disp(sprintf('Calculating features...'));
	disp(sprintf('ChannelValidity: [%d %d %d %d]', MClust_ChannelValidity));
	for iF = 1:length(featureFiles)
		[p n] = fileparts(featureFiles{iF});
		disp(sprintf('..%s', n));
	end
end

Write_fd_file(MClust_FDdn, [MClust_TTdn filesep MClust_TTfn MClust_TText], ...
	featureList, MClust_ChannelValidity, record_block_size, template_matching, NormalizeFDYN)

% Load features 

% Check FD channel validity match
ok = true(size(featureFiles)); temp = cell(size(featureFiles));
for iF = 1:length(featureFiles)
   temp{iF} = load([featureFiles{iF} MClust_FDext],'-mat');
   ok(iF) = all(MClust_ChannelValidity == temp{iF}.ChannelValidity);
end
if any(~ok)
	errordlg({'Channel Validity in feature files'; ...
		      'do not match user-defined Channel Validity.'; ...
			  'Delete problematic feature data files and recalculate them.'},...
		'MClust error', 'modal');
	disp(sprintf('User channel validity: [%d %d %d %d]', MClust_ChannelValidity));
	for iF = 1:length(featureFiles)
		disp(sprintf('%s: [%d %d %d %d]', featureFiles{iF}, temp{iF}.ChannelValidity));
	end
end

% Calculate features into memory

% count features
nFeat = zeros(size(featureList));
nSamps = zeros(size(featureList));
for iF = 1:length(featureList)
	 temp = load([featureFiles{iF} MClust_FDext],'-mat');
	 nFeat(iF) = length(temp.FeatureNames);
	 nSamps(iF) = length(temp.FeatureIndex);
end
if length(unique(nSamps)) > 1
    error('MClust:Error', {'Number of samples in FeatureData files do not match.','Delete feature data files and recompute.'});
end

nFeat = sum(nFeat);

% get times
temp = load([featureFiles{1} MClust_FDext], '-mat', 'FeatureTimestamps');
MClust_FeatureTimestamps = temp.FeatureTimestamps;

% allocate
MClust_FeatureNames = cell(nFeat,1);
MClust_FeatureSources = cell(nFeat,2);
iC = 1; 
for iF = 1:length(featureList)
    fn = [featureFiles{iF} MClust_FDext];
    temp = load(fn,'-mat', 'FeatureNames');
    for iN = 1:length(temp.FeatureNames)
        MClust_FeatureNames{iC} = temp.FeatureNames{iN};
        MClust_FeatureSources{iC,1} = fn;
        MClust_FeatureSources{iC,2} = iN;
        iC = iC+1;
	end
end