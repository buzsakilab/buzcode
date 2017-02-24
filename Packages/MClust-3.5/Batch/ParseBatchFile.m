function [gPar,fPar,def,KeyWord,Values] = ParseBatchFile(fn)
%
% function [gPar,fPar,def,KeyWord,Values] = ParseBatchFile(fn)
%
%  Parse a Batch.txt command file and return parameters in
%  cell arrays gPar (global parametes)  and fPar (file parameters)
%
% Original code by PL 2000.
% Modified by ADR/NCST 2002.
%
% Status: PROMOTED (Release version)
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


AllowedKeyWords = { ...
	'BEGINGLOBALS',...
	'SPIKEFILESSUBDIR',...
	'FINDFILESSEARCHSTRING',...
	'PROCESSINGDIRECTORY',...
	'CLEANEDSPIKEFILESSUBDIR',...
	'CLEANEDSPIKEFILESPREFIX',...
	'CLEANEDSPIKEFILESPOSTFIX',...
	'CLEANEDSPIKEFILESEXT',...
	'FEATUREDATADIR', ... % feature data directory name (FD)
	'LOADINGENGINE', ... % loading engine to use (LoadTT0_NT)
	'USEFEATURES', ... % features to use -- default {'energy','wavePC1','wavePC2','waveFFT'}
	'EXTRAFEATURES', ... % features to pass in but not autocluster on -- default {'peak','valley'}
	'KKWIK_MINCLUSTERS', ... %minimum number of clusters to begin KlustaKwik with
	'KKWIK_MAXCLUSTERS', ... % maximum number of clusters to begin KlustaKwik with
	'SUBSAMPLEDFILESSUBDIR',...
	'SUBSAMPLEDFILESPREFIX',...
	'SUBSAMPLEDFILESPOSTFIX',...
	'SUBSAMPLEDFILESEXT',...
	'USECOINCIDENCEDETECTION',...
	'CLUSTERALGORITHM',...
	'REMOVEFROMBATCHLIST',...
	'ADDTOBATCHLIST',...
	'ENDGLOBALS',...
	'BEGINDEFAULTS',...
	'RECORDINGSYSTEM',...
	'NUMBEROFCHANNELSPERPROBE',...
	'CHANNELVALIDITY',...
	'SUBSAMPLETONSPIKES',...
	'SUBSAMPLEMETHOD',...
	'SUBSAMPLEPARAMETERS',...
	'NUMBEROFNEARESTNEIGHBORS',...
	'ENDDEFAULTS',...
	'BEGINFILES',...
	'FILE',...
	'ENDFILES' ...
	};

%% set defaults
gPar.SearchString = '';
gPar.ProcessingDirectory = '.';
gPar.SpikeFilesSubDir = '.';
gPar.CleanedSpikeFilesSubDir = '';
gPar.CleanedSpikeFilesPrefix = 'TT';
gPar.CleanedSpikeFilesPostfix = '';
gPar.CleanedSpikeFilesExt = '.tt';
gPar.FeatureDataDir = 'FD';
gPar.LoadingEngine = '';
gPar.KKwikMinClusters = 20;
gPar.KKwikMaxClusters = 30;
gPar.UseFeatures = {'energy','wavePC1','wavePC2','waveFFT'};
gPar.ExtraFeatures = {'peak','valley'};
gPar.SubsampledFilesSubDir = '';
gPar.SubsampledFilesPrefix = '';
gPar.SubsampledFilesPostfix = 's';
gPar.SubsampledFilesExt = '.tt';
gPar.RemoveList = {''};
gPar.AddList = {};
gPar.FileList = {};
gPar.CleanedFileNames = {};
gPar.SubsampledFileNames = {};
gPar.FeatureDataFileNames = {};
gPar.UseCoincidenceDetection = 0;
gPar.ClusterAlgorithm = 'BBClust';

def.FileName = '';
def.RecordingSystem = 'NT';
def.NumberOfChannels = 4;
def.ChannelValidity = [1 1 1 1];
def.SubsampleToNSpikes = 0;
def.SubsampleMethod = 1;
def.SubsampleParameters = 0.8;
def.NN = 20;
def.MinClusters = 20;
def.MaxClusters = 30;

fPar = {};

if nargin == 0
	fn = 'Batch.txt';
end

% load Batch.txt into a cellarray with one word or parameter per cell
if version('-release') >= 12
	lines = textread(fn,'%s','commentstyle','matlab');
else
	lines = textread(fn,'%s','delimiter','\n','commentstyle','matlab');
end%if
blanklines(strcmp(lines, '')) = [];  % remove cells with empty strings

% we now have a list of keywords and the parameters following each keyword (if any)
% convert into two parallel cell arrays KeyWord{


iKeyWord = 0;
KeyWord = {};
Values = {};
IsAKeyWord = 0;
for ii=1:length(lines)
	% decide if line is a Keyword or a Value depending on its '--' characters
	if length(lines{ii}) < 2
		IsAKeyWord = 0;              % line is certainly not a keyword
	else
		if strcmp(lines{ii}(1:2),'--')  % line is definitely a keyword
			IsAKeyWord = 1;
		else
			IsAKeyWord = 0;
		end
	end
	% fill parallel cell arrays KeyWord and Values
	if IsAKeyWord
		% lines(ii) is a keyword
		iKeyWord = iKeyWord + 1;
		KeyWord{iKeyWord} = lines{ii}(3:end);         % strip off the beginning -- from each keyword
		Values{iKeyWord} = {};
		iValue = 0;
	else
		% lines(ii) is a value to the current keyword
		iValue = iValue + 1;
		Values{iKeyWord}{iValue} = lines{ii};
	end
end

nKeyWords = length(KeyWord);


%% process keywords
KeyWord = upper(KeyWord);    %% convert keywords to upper

%check if all keywords are in the list of reckognized keywords
for ii = 1:nKeyWords
	isAmatch = strmatch(KeyWord{ii},AllowedKeyWords,'exact');
	if isempty(isAmatch)
		error( ['Keyword ' KeyWord{ii} ' is not a reckognized KeyWord.' ...
			' Check correct spelling or inclusion of' ...
			' an accidental blank between the doubledash -- and a keyword!']);
	end
end

%give all empty Values a default empty string
for ii = 1:nKeyWords
	if isempty(Values{ii})
		Values{ii}{1} = '';
	end
end



%process global parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find index range for global paramters
ig_first = strmatch('BEGINGLOBALS',KeyWord,'exact');
ig_last  = strmatch('ENDGLOBALS',KeyWord,'exact');
if isempty(ig_first) || isempty(ig_last)
	error('MClust:Batch', 'Could not find the BEGINGLOBALS or ENDGLOBALS keyword. These are required exactly once!');
end
if length(ig_first)>1 || length(ig_last)>1
	error('MClust:Batch', 'There are multiple BEGINGLOBALS or ENDGLOBALS keywords. These are required exactly once!');
end

%extract keywords and values
for ii = ig_first+1 : ig_last-1

	switch KeyWord{ii}

		case 'FINDFILESSEARCHSTRING'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.SearchString = Values{ii}{1};

		case 'SPIKEFILESSUBDIR'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.SpikeFilesSubDir = Values{ii}{1};

		case 'PROCESSINGDIRECTORY'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.ProcessingDirectory= Values{ii}{1};

		case 'CLEANEDSPIKEFILESSUBDIR'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.CleanedSpikeFilesSubDir = Values{ii}{1};

		case 'CLEANEDSPIKEFILESPREFIX'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.CleanedSpikeFilesPrefix = Values{ii}{1};

		case 'CLEANEDSPIKEFILESPOSTFIX'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.CleanedSpikeFilesPostfix = Values{ii}{1};

		case 'CLEANEDSPIKEFILESEXT'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.CleanedSpikeFilesExt = Values{ii}{1};

		case 'SUBSAMPLEDFILESSUBDIR'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.SubsampledFilesSubDir = Values{ii}{1};

		case 'FEATUREDATADIR'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.FeatureDataDir = Values{ii}{1};

		case 'LOADINGENGINE'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.LoadingEngine = Values{ii}{1};

		case 'USEFEATURES'
			if isempty(Values{ii})
				error( ['Keyword ' KeyWord{ii} ' takes one or more parameter!'] );
			end
			gPar.UseFeatures = Values{ii};

		case 'EXTRAFEATURES'
			if isempty(Values{ii})
				error( ['Keyword ' KeyWord{ii} ' takes one or more parameter!'] );
			end
			gPar.ExtraFeatures = Values{ii};

		case 'SUBSAMPLEDFILESPREFIX'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.SubsampledFilesPrefix = Values{ii}{1};

		case 'SUBSAMPLEDFILESPOSTFIX'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.SubsampledFilesPostfix = Values{ii}{1};

		case 'SUBSAMPLEDFILESEXT'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			gPar.SubsampledFilesExt = Values{ii}{1};

		case 'USECOINCIDENCEDETECTION'
			if length(Values{ii}) ~= 1
				error( ['Keyword ' KeyWord{ii} ' takes exactly one parameter 0 or 1!'] );
			end
			gPar.UseCoincidenceDetection = str2num(Values{ii}{1});

		case 'CLUSTERALGORITHM'
			if length(Values{ii}) ~= 1
				error( ['Keyword ' KeyWord{ii} ' takes exactly one parameter 0 or 1!'] );
			end
			switch Values{ii}{1}
				case 'BBClust'
					error('MClust:Batch', 'BBclust is no longer supported in MClust.');
				case 'KlustaKwik'
					gPar.ClusterAlgorithm = Values{ii}{1};
				otherwise
					error('MClust:Batch','Only KlustaKwik is currently supported in MClust Batch processing.');
			end%switch

		case 'REMOVEFROMBATCHLIST'
			gPar.RemoveList = [gPar.RemoveList, Values{ii}];   %append to previous list

		case 'ADDTOBATCHLIST'
			gPar.AddList = [gPar.AddList, Values{ii}];   %append to previous list

	end%switch

end%for ii

% goto ProcessingDirectory
pushdir(gPar.ProcessingDirectory);

% construct gPar.FileList
fnames = sort(gPar.AddList)';
fnames_tmp = {};
if ~isempty(gPar.SearchString)
	fnames_tmp = FindFiles(gPar.SearchString,...
		'StartingDirectory',gPar.SpikeFilesSubDir, ...
		'CheckSubdirs',0);
	%fnames_tmp = FindFiles(gPar.SearchString); % ncst debug 12/05/01
	for ii = 1:length(fnames_tmp)
		[path, fnames_tmp{ii}, ext] = fileparts(fnames_tmp{ii});     %strip off dirpath
		fnames_tmp{ii} = [fnames_tmp{ii} ext];                       %add extension
	end
	fnames = [fnames; fnames_tmp];
end

jj = 0;
fnames_tmp = {};
for ii = 1:length(fnames)
	if exist([gPar.SpikeFilesSubDir filesep fnames{ii}],'file') == 0
		error([' ParseBatchFile:  file ', fnames{ii}, ' does not exist! Typo???' ]);
	else
		jj = jj+1;
		[path, fnames_tmp{jj}, ext] = fileparts(fnames{ii});     %strip off dirpath
		fnames_tmp{jj} = [fnames_tmp{jj} ext];                       %add extension
	end%if
end


gPar.FileList = unique(fnames_tmp);
if isempty(gPar.FileList)
	error('FileList is empty: No spike files to process');
end%if
for ii = 1:length(gPar.RemoveList)
	indx = strmatch(lower(gPar.RemoveList{ii}), lower(gPar.FileList),'exact');
	if ~isempty(indx)
		gPar.FileList(indx) = [];
	elseif ~isempty(gPar.RemoveList{ii})
		warning('MClust:Batch', ['Could not remove file ' gPar.RemoveList{ii} ' from file list. File not found in list!']);
	end
end%for



%process default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find index range for global paramters
ig_first = strmatch('BEGINDEFAULTS',KeyWord,'exact');
ig_last  = strmatch('ENDDEFAULTS',KeyWord,'exact');
if isempty(ig_first) || isempty(ig_last)
	error('MClust:Batch', 'Could not find the BEGINDEFAULTS or ENDDEFAULTS keyword. These are required exactly once!');
end
if length(ig_first)>1 || length(ig_last)>1
	error('MClust:Batch','There are multiple BEGINDEFAULTS or ENDDEFAULTS keywords. These are required exactly once!');
end

%extract keywords and values
for ii = ig_first+1 : ig_last-1

	switch KeyWord{ii}

		case 'RECORDINGSYSTEM'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			end
			vv = Values{ii}{1};
			if strcmpi(vv,'NT')
				def.RecordingSystem = 'NT';
			elseif strcmpi(vv,'SUN')
				def.RecordingSystem = 'SUN';
			else
				error(['Only SUN and NT are allowed as values for ' KeyWord{ii}]);
			end%if

		case 'NUMBEROFCHANNELSPERPROBE'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			elseif length(Values{ii}) == 1
				def.NumberOfChannels = str2num(Values{ii}{1});
			end

		case 'CHANNELVALIDITY'
			if length(Values{ii}) ~= 4
				error( ['Keyword ' KeyWord{ii} ' takes exactly 4 parameters!'] );
			end
			for jj=1:4
				def.ChannelValidity(jj) = str2num(Values{ii}{jj});
			end

		case 'SUBSAMPLETONSPIKES'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			elseif length(Values{ii}) == 1
				def.SubsampleToNSpikes = str2num(Values{ii}{1});
			end

		case 'SUBSAMPLEMETHOD'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			elseif length(Values{ii}) == 1
				def.SubsampleMethod = str2num(Values{ii}{1});
			end

		case 'SUBSAMPLEPARAMETERS'
			if length(Values{ii}) > 2
				error( ['Keyword ' KeyWord{ii} ' takes 0,1 or 2 parameters!'] );
			elseif ~isempty(Values{ii})
				for jj=1:length(Values{ii})
					def.SubsampleParameters(jj) = str2num(Values{ii}{jj});
				end
			end

		case 'KKWIK_MINCLUSTERS'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			elseif length(Values{ii}) == 1
				def.KKwikMinClusters = str2num(Values{ii}{1});
			end

		case 'KKWIK_MAXCLUSTERS'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			elseif length(Values{ii}) == 1
				def.KKwikMaxClusters = str2num(Values{ii}{1});
			end

		case 'NUMBEROFNEARESTNEIGHBORS'
			if length(Values{ii}) > 1
				error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
			elseif length(Values{ii}) == 1
				def.NN = str2num(Values{ii}{1});
			end

	end%switch
end%for


%process individual file parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set defaults
for ii = 1:length(gPar.FileList)
	fPar{ii} = def;
	fPar{ii}.FileName = gPar.FileList{ii};
end%for

%find indexes of --FILE keywords
ix = strmatch('FILE',KeyWord,'exact');
if isempty(ix)
	%warning( ['Could not find any FILE keyword. Only defaults apply!']);
end

%loop over all --FILE keywords
for kk = 1:length(ix)
	ig_first = ix(kk)+1;
	if kk < length(ix)
		ig_last = ix(kk+1)-1;
	else
		ig_last = length(KeyWord);
	end

	% loop over all files in current --FILE section
	flist = Values{ix(kk)};
	for jj = 1:length(flist)
		%find index for current file in gPar.FileList
		findx = strmatch(lower(flist(jj)),lower(gPar.FileList),'exact');
		if isempty(findx)
			error(['File ' flist{jj} ' does not exist! Remove it from --File list!']);
		end
		if length(indx) > 1
			error(['Duplicate file name ' flist{jj} ' in File List!']);
		end

		%extract keywords and values and override default for current file
		for ii = ig_first : ig_last

			switch KeyWord{ii}

				case 'RECORDINGSYSTEM'
					if length(Values{ii}) > 1
						error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
					end
					vv = Values{ii}{1};
					if strcmpi(vv,'NT')
						fPar{findx}.RecordingSystem = 'NT';
					elseif strcmpi(vv,'SUN')
						fPar{findx}.RecordingSystem = 'SUN';
					else
						error(['Only SUN and NT are allowed as values for ' KeyWord{ii}]);
					end%if

				case 'KKWIK_MINCLUSTERS'
					if length(Values{ii}) > 1
						error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
					elseif length(Values{ii}) == 1
						fPar{findx}.KKwikMinClusters = str2num(Values{ii}{1});
					end

				case 'KKWIK_MAXCLUSTERS'
					if length(Values{ii}) > 1
						error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
					elseif length(Values{ii}) == 1
						fPar{findx}.KKwikMaxClusters = str2num(Values{ii}{1});
					end

				case 'NUMBEROFCHANNELSPERPROBE'
					if length(Values{ii}) > 1
						error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
					elseif length(Values{ii}) == 1
						fPar{findx}.NumberOfChannels = str2num(Values{ii}{1});
					end

				case 'CHANNELVALIDITY'
					if length(Values{ii}) ~= 4
						error( ['Keyword ' KeyWord{ii} ' takes exactly 4 parameters!'] );
					end
					for iCV=1:4
						fPar{findx}.ChannelValidity(iCV) = str2num(Values{ii}{iCV});
					end

				case 'SUBSAMPLETONSPIKES'
					if length(Values{ii}) > 1
						error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
					elseif length(Values{ii}) == 1
						fPar{findx}.SubsampleToNSpikes = str2num(Values{ii}{1});
					end

				case 'SUBSAMPLEMETHOD'
					if length(Values{ii}) > 1
						error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
					elseif length(Values{ii}) == 1
						fPar{findx}.SubsampleMethod = str2num(Values{ii}{1});
					end

				case 'SUBSAMPLEPARAMETERS'
					if length(Values{ii}) > 2
						error( ['Keyword ' KeyWord{ii} ' takes 0,1 or 2 parameters!'] );
					elseif ~isempty(Values{ii})
						for iSP=1:length(Values{ii})
							fPar{findx}.SubsampleParameters(iSP) = str2num(Values{ii}{iSP});
						end
					end

				case 'NUMBEROFNEARESTNEIGHBORS'
					if length(Values{ii}) > 1
						error( ['Keyword ' KeyWord{ii} ' takes none or one parameter!'] );
					elseif length(Values{ii}) == 1
						fPar{findx}.NN = str2num(Values{ii}{1});
					end

			end%switch
		end%for ii

	end%for jj

end%for kk


% go back to previous directory
popdir;
