function [lfp] = bz_GetLFP(varargin)
% bz_GetLFP - Get local field potentials.
%
%  Load local field potentials from disk. No longer dependent on
%  FMAT/SetCurrentSession.
%
%  USAGE
%
%    [lfp] = bz_GetLFP(channels,<options>)
%
%  INPUTS
%
%    channels(required) -must be first input, numeric  
%                        list of channels to load (use keyword 'all' for all)
%                        channID is 0-indexing, a la neuroscope
%  Name-value paired inputs:
%    'basepath'           - folder in which .lfp file will be found (default
%                           is pwd)
%                           folder should follow buzcode standard:
%                           whateverPath/baseName
%                           and contain file baseName.lfp
%    'basename'           -base file name to load
%    'intervals'          -list of time intervals [0 10; 20 30] to read from 
%                           the LFP file (default is [0 inf])
%    'downsample'         -factor to downsample the LFP (i.e. 'downsample',5
%                           will load a 1250Hz .lfp file at 250Hz)
%    'noPrompts'          -logical (default) to supress any user prompts
%
%  OUTPUT
%
%    lfp             struct of lfp data. Can be a single struct or an array
%                    of structs for different intervals.  lfp(1), lfp(2),
%                    etc for intervals(1,:), intervals(2,:), etc
%    .data           [Nt x Nd] matrix of the LFP data
%    .timestamps     [Nt x 1] vector of timestamps to match LFP data
%    .interval       [1 x 2] vector of start/stop times of LFP interval
%    .channels       [Nd X 1] vector of channel ID's
%    .samplingRate   LFP sampling rate [default = 1250]
%    .duration       duration, in seconds, of LFP interval
%
%
%  EXAMPLES
%
%    % channel ID 5 (= # 6), from 0 to 120 seconds
%    lfp = bz_GetLFP(5,'restrict',[0 120]);
%    % same, plus from 240.2 to 265.23 seconds
%    lfp = bz_GetLFP(5,'restrict',[0 120;240.2 265.23]);
%    % multiple channels
%    lfp = bz_GetLFP([1 2 3 4 10 17],'restrict',[0 120]);
%    % channel # 3 (= ID 2), from 0 to 120 seconds
%    lfp = bz_GetLFP(3,'restrict',[0 120],'select','number');

% Copyright (C) 2004-2011 by Michaël Zugaro
% editied by David Tingley, 2017
%
% added support for NWB by Konstantinos Nasiotis, 2019

% NOTES
% -'select' option has been removed, it allowed switching between 0 and 1
%   indexing.  This should no longer be necessary with .lfp.mat structs
%
% TODO
% add saveMat input 
% expand channel selection options (i.e. region or spikegroup)
% add forcereload
%% Parse the inputs!

channelsValidation = @(x) isnumeric(x) || strcmp(x,'all');

% parse args
p = inputParser;
addRequired(p,'channels',channelsValidation)
addParameter(p,'basename','',@isstr)
addParameter(p,'intervals',[],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'downsample',1,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);

parse(p,varargin{:})
basename = p.Results.basename;
channels = p.Results.channels;
downsamplefactor = p.Results.downsample;
basepath = p.Results.basepath;
noPrompts = p.Results.noPrompts;

% doing this so you can use either 'intervals' or 'restrict' as parameters to do the same thing
intervals = p.Results.intervals;
restrict = p.Results.restrict;
if isempty(intervals) && isempty(restrict) % both empty
    intervals = [0 Inf];
elseif isempty(intervals) && ~isempty(restrict) % intervals empty, restrict isn't
    intervals = restrict;
end

%% let's check that there is an appropriate LFP file
if isempty(basename)
   %disp('No basename given, so we look for a *lfp/*eeg file...')
   d = dir([basepath filesep '*lfp']);
   if length(d) > 1 % we assume one .lfp file or this should break
       error('there is more than one .lfp file in this directory?');
   elseif length(d) == 0
       d = dir([basepath filesep '*eeg']);
       if isempty(d)
           % If no .lfp file is present, check for a .nwb
           d = dir([basepath filesep '*.nwb']);
           if length(d) > 1 % If more than one .nwb files exist in the directory, select which one to load from
               warning('there is more than one .nwb file in this directory');
               % Check if a specific behavior was called to be loaded. If not, display a
               % pop-up list with the available behaviors for selection
               [iNWBFile, ~] = listdlg('PromptString','Which NWB file would you like to load?',...
                                         'ListString',{d.name},'SelectionMode','single');
               d = d(iNWBFile);
           else
               if isempty(d)
                    error('could not find an lfp/eeg/nwb file..')
               end
           end
       end
   end
   lfp.Filename = d.name;
   basename = strsplit(lfp.Filename,'.');
   if length(basename) > 2
       base = [];
       for i=1:length(basename)-1
          base = [base basename{i} '.'];
       end
       basename = base(1:end-1);  % this is an fugly hack to make things work with Kenji's naming system...
   else
       basename = basename{1};
   end
   
else
   d = dir([basepath filesep basename '.lfp']);
   if length(d) > 1 % we assume one .lfp file or this should break
       error('there is more than one .lfp file in this directory?');
   elseif length(d) == 0
       d = dir([basepath filesep basename '.eeg']);
       if isempty(d)
           d = dir([basepath filesep basename '.nwb']);
           if length(d) > 1 % we assume one .nwb file or this should break
               warning('there is more than one .nwb file in this directory');
               % Check if a specific behavior was called to be loaded. If not, display a
               % pop-up list with the available behaviors for selection
               [iNWBFile, ~] = listdlg('PromptString','Which NWB file would you like to load?',...
                                         'ListString',{d.name},'SelectionMode','single');
               d = d(iNWBFile);
           else
               if isempty(d)
                    error('could not find an lfp/eeg/nwb file..')
               end
           end
       end
   end
   lfp.Filename = d.name;   
end


%% Add a flag that shows if variables are loaded from a .nwb file
%  nwb has everything saved in a single file

if ~isempty(strfind(lfp.Filename,'nwb'))
    nwb_loaded = 1;
else
    nwb_loaded = 0;
end

%% things we can parse from sessionInfo or xml file

if ~nwb_loaded
    sessionInfo = bz_getSessionInfo(basepath, 'noPrompts', noPrompts);
else
    % Load the .nwb info
    nwb2 = nwbRead(lfp.Filename);
    sessionInfo = get_SessionInfoFromNWB(nwb2); % The get_SessionInfoFromNWB loads only the needed fields for the functions tested. Maybe improve
    
    if ~sessionInfo.LFPDataPresent
        error('No LFP data detected in this file. Make sure that there is a field within: nwb2.processing.get("ecephys").nwbdatainterface.get("LFP").electricalseries')
    end
end

try
    samplingRate = sessionInfo.lfpSampleRate;
catch
    samplingRate = sessionInfo.rates.lfp; % old ugliness we need to get rid of
end
samplingRateLFP = samplingRate./downsamplefactor;

if mod(samplingRateLFP,1)~=0
    error('samplingRate/downsamplefactor must be an integer')
end
%% Channel load options
%Right now this assumes that all means channels 0:nunchannels-1 (neuroscope
%indexing), we could also add options for this to be select region or spike
%group from the xml...
if strcmp(channels,'all')
    channels = sessionInfo.channels;
else
    %display('Loading Channels ',num2str(channels),' (0-indexing, a la Neuroscope)')
end

%% get the data
disp('loading LFP file...')
nIntervals = size(intervals,1);
% returns lfp/bz format
for i = 1:nIntervals
    lfp(i).duration = (intervals(i,2)-intervals(i,1));
    lfp(i).interval = [intervals(i,1) intervals(i,2)];

    % Load data and put into struct
    % we assume 0-indexing like neuroscope, but bz_LoadBinary uses 1-indexing to
    % load....
    
    if ~nwb_loaded
        lfp(i).data = bz_LoadBinary([basepath filesep lfp.Filename],...
            'duration',double(lfp(i).duration),...
                      'frequency',samplingRate,'nchannels',sessionInfo.nChannels,...
                      'start',double(lfp(i).interval(1)),'channels',channels+1,...
                      'downsample',downsamplefactor);
                  
    else
        %% NWB
        % Memory preallocation is used for faster assignment
        if intervals(i,2) == Inf
            use_this_bracket = sessionInfo.samples_NWB/sessionInfo.lfpSampleRate;
            data_temp = zeros(length(channels),ceil((use_this_bracket - intervals(i,1))*samplingRateLFP));
        else
            use_this_bracket = intervals(i,2);
            data_temp = zeros(length(channels),floor((intervals(i,2) - intervals(i,1))*samplingRateLFP));
        end

        % Differentiate between sequential loading of each channel and
        % simultaneous loading of neighboring channels
        if all(diff(channels) == 1)
        	data_temp = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(sessionInfo.LFPDataKey).data.load([channels(1)+1, intervals(i,1)*samplingRate+1],[1, downsamplefactor],[channels(end)+1, use_this_bracket*samplingRate]);
        else
            for iChannel = 1:length(channels)
                disp(['Now concatenating channel: ' num2str(channels(iChannel))])
                data_temp(iChannel,:) = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(sessionInfo.LFPDataKey).data.load([channels(iChannel)+1, intervals(i,1)*samplingRate+1],[1, downsamplefactor],[channels(iChannel)+1, use_this_bracket*samplingRate]);
            end 
        end
        lfp(i).data = data_temp'; clear data_temp
    end
    
    lfp(i).timestamps = [lfp(i).interval(1):(1/samplingRateLFP):...
                        (lfp(i).interval(1)+(length(lfp(i).data)-1)/...
                        samplingRateLFP)]';
    lfp(i).channels = channels;
    lfp(i).samplingRate = samplingRateLFP;
    % check if duration is inf, and reset to actual duration...
    if lfp(i).interval(2) == inf
        lfp(i).interval(2) = length(lfp(i).timestamps)/lfp(i).samplingRate;
        lfp(i).duration = (lfp(i).interval(i,2)-lfp(i).interval(i,1));
    end
    
    if isfield(sessionInfo,'region') && isfield(sessionInfo,'channels')
        [~,~,regionidx] = intersect(lfp(i).channels,sessionInfo.channels,'stable');
        lfp(i).region = sessionInfo.region(regionidx); % match region order to channel order..
    end
end
