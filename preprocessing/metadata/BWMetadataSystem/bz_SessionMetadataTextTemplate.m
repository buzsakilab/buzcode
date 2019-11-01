function SessionMetadata = bz_SessionMetadataTextTemplate(basepath)
% This .m file is a generic template and will be copied to session
% folders/session folders and will be called basename_SessionMetadataText.m.
% It is represents code that will be used to
% create a basename.SessionMetadata.mat file.  Edits here will be reflected
% in the .SessionMetadata.mat created after this file is run.  
% 
% Developers interested in changing how SessionMetadata is made generally 
% should edit bz_SessionMetadataTextTemplate.m 
%
% INPUTS
%   basepath - computer path to the session folder of interest
%   AnimalSessionUncertain - If 1, user must confirm source of
%      animalmetadata.mat. If 0 will look for it.  Default = 1
% 
% OUTPUTS
%    SessionMetadata and saved basename.SessionMetadata.mat
%       - Struct array containing fields including number of drive turns
%       since surgery, bad channels, behavioral info etc.
%
% See also: bz_EditSessionMetadata, bz_RunSessionMetadata, bz_AnimalMetadataTextTemplate
%
% Brendon Watson 2017


%% Input handling/parsing
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);


%% HUMAN INPUT THIS SECTION
% Preprocessing specifiers
SessionMetadata.Preprocess.AnimalMetadataSource = {'AnimalFolderAbove','Basepath'};%!blank will lead to user being asked to find file

% Surgery and Animal metadata
SessionMetadata.Animal.WeightGrams = [];

% Extracell Ephys metadata
SessionMetadata.ExtracellEphys.NumberOfTurnsSinceSurgery = [0];%vector, one entry per probe
SessionMetadata.ExtracellEphys.Probes.PluggingOrder = ['FromAnimalMetadata'];%vector, one entry per probe. blank defaults to animal plugging order

% SessionMetadata.ExtracellEphys.BadShanks = [];% vector for this recording. base 1.  NOT REALLY USED YET
     % These bad shanks will be used to populate bad channels
SessionMetadata.ExtracellEphys.BadChannels = 'FromAnimalMetadata';%... or number base 0
SessionMetadata.ExtracellEphys.ChannelNotes = {''};


%if these values are "FromAnimalMetadata", they'll be auto-populated w
%numbers later below. Otherwise a number entered here will be used
% SessionMetadata.ExtracellEphys.NumberOfChannels = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.SpikeGroups = {'FromAnimalMetaData','FromXML','FromSessionInfo'}; %Pick one, delete the other.

%For setting up to record extra channels for opto, frames etc.  
% NumExtraChansPerExtraGroup field will be a vector of number of channels per
% extra group.   The number of groups is set by the number of entries... 
% the number within each group is signified by the number in each
% entry.  (ie [3 5] signifies 2 post-probe groups, first with 3 channels, the 
% second with 5 channels.
SessionMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.ExtraChannels.GroupNames = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.ExtraChannels.ChannelNames = 'FromAnimalMetadata';%... or number

%Basic Ephys params
SessionMetadata.ExtracellEphys.Parameters.SampleRate = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.Parameters.Amplification = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.Parameters.VoltsPerUnit = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.Parameters.BitsPerSample = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.Parameters.VoltageRange = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.Parameters.PointsPerWaveform = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.Parameters.PeakPointInWaveform = 'FromAnimalMetadata';%... or number
SessionMetadata.ExtracellEphys.Parameters.FeaturesPerWave = 'FromAnimalMetadata';%... or number



% Optogenetics metadata
% FiberNum = 1;%copy more of these blocks - one per probe
% SessionMetadata.Optogenetics.Fibers(FiberNum).StimulusRecordingChannelDatType = 'AnalogIn';%'AnalogIn' for Intan, 'Main' if fused with other dat
% SessionMetadata.Optogenetics.Fibers(FiberNum).StimulusRecordingChannelNumber = [];%
%?else

% Behavior metadata - this section makes an error
% BehaviorSessionNumber = 1;
% SessionMetadata.BehaviorEpisodes.(BehaviorSessionNumber).StartStopSeconds = [];%startstop pair, in seconds
% SessionMetadata.BehaviorEpisodes(BehaviorSessionNumber).DatFileIndex = [];%which dat files are included in this Episode
% SessionMetadata.BehaviorEpisodes(BehaviorSessionNumber).BehaviorType = '';%verbal description of behavior type, ie "LinearTrack"
% SessionMetadata.BehaviorEpisodes(BehaviorSessionNumber).BehaviorNotes = '';
% SessionMetadata.BehaviorEpisodes.(BehaviorSessionNumber).Motion.MotionTrackerType = '';%Optitrack or Basler etc
% SessionMetadata.BehaviorEpisodes.(BehaviorSessionNumber).Motion.CameraFrameChannelType = 'DigitalIn';%'DigitalIn' for Intan, 'Main' if fused with other dat
% SessionMetadata.BehaviorEpisodes.(BehaviorSessionNumber).Motion.CameraFrameChannelNumber = [];%channel number where frame times are recorded
% SessionMetadata.BehaviorEpisodes(BehaviorSessionNumber).Location = '';
% 
% SleepSessionNumber = 1;
% SessionMetadata.SleepEpisodes.(SleepSessionNumber).StartStopSeconds = [];%startstop pair, in seconds
% SessionMetadata.SleepEpisodes(SleepSessionNumber).DatFileIndex = [];%which dat files are included in this Episode
% SessionMetadata.SleepEpisodes(SleepSessionNumber).BehaviorNotes = '';
% SessionMetadata.SleepEpisodes.(SleepSessionNumber).Motion.MotionTrackerType = '';%Optitrack or Basler etc
% SessionMetadata.SleepEpisodes.(SleepSessionNumber).Motion.CameraFrameChannelType = 'DigitalIn';%'DigitalIn' for Intan, 'Main' if fused with other dat
% SessionMetadata.SleepEpisodes.(SleepSessionNumber).Motion.CameraFrameChannelNumber = [];%channel number where frame times are recorded
% SessionMetadata.SleepEpisodes(SleepSessionNumber).Location = '';

% Virus metadata
    %need to fill this in

% Intracell ephys metadata
    %need to fill this in

% Imaging metadata
    %need to fill this in

% Other for ?
% SessionMetadata.Other = {''};


%% Automated after this point, depending on modules used
%% Find AnimalMetadata in order to open it.  
% copy into this folder - this will give crucial info for subsequent work here
% Make an initial version of output from this, overwrite with actual
% metadata from recording system
% 
supradir = fileparts(basepath);
finalamdpath = fullfile(basepath,[basename,'.AnimalMetadata.mat']);
% if we know the AnimalMetadata is one folder above basepath (ie in Animal
% folder as usual)
if strcmp(lower(SessionMetadata.Preprocess.AnimalMetadataSource{1}),'animalfolderabove')
    d2 = dir(fullfile(supradir,'*.AnimalMetadata.mat'));
    inputamdpath = fullfile(supradir,d2(1).name);
elseif strcmp(lower(SessionMetadata.Preprocess.AnimalMetadataSource{1}),'basepath')
    d2 = dir(fullfile(basepath,'*.AnimalMetadata.mat'));
    inputamdpath = fullfile(basepath,d2(1).name);
end

% in other cases search for it and ask user
if ~exist('inputamdpath','var')
    guessamdpath = [];
    d1 = dir(fullfile(basepath,'*.AnimalMetadata.mat'));
    d2 = dir(fullfile(supradir,'*.AnimalMetadata.mat'));
    finalamdpath = fullfile(basepath,[basename,'.AnimalMetadata.mat']);


    if ~isempty(d1)%if already one in this path
        guessamdpath = fullfile(basepath,d1(1).name);
    elseif ~isempty(d2)%if one in directory above
        guessamdpath = fullfile(supradir,d2(1).name);
    else%if cannot be found, ask user 
        disp('No local AnimalMetadata.mat found');
    end
    inputamdpath = fullfile(PathName,FileName);
end

if ~strcmp(inputamdpath,finalamdpath)
    copyfile(inputamdpath,finalamdpath)
end

load(finalamdpath);
SessionMetadata.AnimalMetadata = AnimalMetadata;

clear d1 d2 supradir upperamdpath FileName PathName inputamdpath finalamdpath


%%
if SessionMetadata.AnimalMetadata.Modules.ExtracellEphys
    SessionMetadata.ExtracellEphys.Probes.ProbeDepths = SessionMetadata.ExtracellEphys.NumberOfTurnsSinceSurgery .* SessionMetadata.AnimalMetadata.ExtracellEphys.Probes.UmPerScrewTurn;
    SessionMetadata.ExtracellEphys.Probes.ProbeCenterCoordinates = [];%someone should do this
    SessionMetadata.ExtracellEphys.Probes.ChannelCenterCoordinates = [];%someone should do this
    
    % Read recording system-based metadata from either .rhd or .meta (Intan or Amplirec)
    % Check for compatibility with AnimalMetadata, overwrite using this and
    % warn the user of the conflict.
    SessionMetadata.ExtracellEphys.DatTimingInfo = bz_DatFileMetadata(basepath);

    %now things regardless of input system

    %take from animal metadata if user said so (default for now 3/2019)
    if strcmp(lower(SessionMetadata.ExtracellEphys.BadChannels),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.BadChannels = AnimalMetadata.ExtracellEphys.Channels.BadChannels;
    end
    
    if strcmp(lower(SessionMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup = AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup;
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.ExtraChannels.GroupNames),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.ExtraChannels.GroupNames = AnimalMetadata.ExtracellEphys.ExtraChannels.GroupNames;
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.ExtraChannels.ChannelNames),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.ExtraChannels.ChannelNames = AnimalMetadata.ExtracellEphys.ExtraChannels.ChannelNames;
    end
    
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.SampleRate),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.SampleRate = AnimalMetadata.ExtracellEphys.Parameters.SampleRate;%
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.Amplification),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.Amplification = AnimalMetadata.ExtracellEphys.Parameters.Amplification;%
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.VoltsPerUnit),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.VoltsPerUnit = AnimalMetadata.ExtracellEphys.Parameters.VoltsPerUnit;%
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.BitsPerSample),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.BitsPerSample = AnimalMetadata.ExtracellEphys.Parameters.BitsPerSample;%
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.VoltageRange),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.VoltageRange = AnimalMetadata.ExtracellEphys.Parameters.VoltageRange;%
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate = AnimalMetadata.ExtracellEphys.Parameters.LfpSampleRate;%
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.PointsPerWaveform),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.PointsPerWaveform = AnimalMetadata.ExtracellEphys.Parameters.PointsPerWaveform;%
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.PeakPointInWaveform),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.PeakPointInWaveform = AnimalMetadata.ExtracellEphys.Parameters.PeakPointInWaveform;%
    end
    if strcmp(lower(SessionMetadata.ExtracellEphys.Parameters.FeaturesPerWave),'fromanimalmetadata')
        SessionMetadata.ExtracellEphys.Parameters.FeaturesPerWave = AnimalMetadata.ExtracellEphys.Parameters.FeaturesPerWave;%
    end
    
    %Check for conflicts between original animal-based ephys params and
    %those based on files in this particular recording.
    aparams = SessionMetadata.AnimalMetadata.ExtracellEphys.Parameters;
    sparams = SessionMetadata.ExtracellEphys.Parameters;
    fn = fieldnames(aparams);
    divergentfields = {};
    for fidx = 1:length(fn)
        if getfield(aparams,fn{fidx}) ~= getfield(sparams,fn{fidx})
            divergentfields{end+1} = fn{fidx};
        end
    end
    if ~isempty(divergentfields) % if conflicts, warn user
        warning('Conflicts between ephys parameters specified for Animal and those found in this recording')
        warning('Conflicts found in:')
        disp(divergentfields)
    end
    SessionMetadata.ExtracellEphys.ParametersDivergentFromAnimalMetadata = divergentfields;
    
    
    
%     %Get badchannels from badshanks... ie all channels on bad shanks are
%     %bad
%     % Add bad channels to channels derived from bad channels
% %     if isempty(SessionMetadata.ExtracellEphys.BadChannels)
%         if ~isempty(SessionMetadata.ExtracellEphys.BadShanks)
%             chanpergp = SessionMetadata.AnimalMetadata.ExtracellEphys.Probes.ProbeSpikeGroupLayoutSuperficialToDeep;
%             badchans = [];
%             for sidx = 1:length(SessionMetadata.ExtracellEphys.BadShanks)
%                 tshank = SessionMetadata.ExtracellEphys.BadShanks;
%                 badchans = cat(1,badchans,chanpergp{tshank});
%             end
%             SessionMetadata.ExtracellEphys.BadChannels = cat(1,SessionMetadata.ExtracellEphys.BadChannels(:),badchans(:));
%         end
% %     end 
    
    %Get Spikegroups
    if strcmp(SessionMetadata.ExtracellEphys.SpikeGroups,'FromXML')
        display('Getting spikegroups from the .xml...')
        %Load the XML for spike groups
        xmlname = fullfile(basepath,[basename,'.xml']);
        if ~exist(xmlname,'file')
            [FileName,PathName] = uigetfile('.xml','Find the .xml to load SpikeGroups from',basepath);
            xmlname = fullfile(PathName,FileName);
        end
        xmlparms = LoadParameters(xmlname);
        SessionMetadata.ExtracellEphys.SpikeGroups = xmlparms.SpkGrps;
    elseif strcmp(SessionMetadata.ExtracellEphys.SpikeGroups,'FromSessionInfo')
        if exist(fullfile(basepath,[basename '.SessionInfo.mat']))
            si = bz_getSessionInfo(basepath);
            SessionMetadata.ExtracellEphys.SpikeGroups = si.spikeGroups;
        else
            disp('No SessionInfo found in folder.  Taking spike groups from AnimalMetadata')
            SessionMetadata.ExtracellEphys.SpikeGroups = AnimalMetadata.ExtracellEphys.SpikeGroups;
        end
    elseif strcmp(SessionMetadata.ExtracellEphys.SpikeGroups,'FromAnimalMetaData')
        SessionMetadata.ExtracellEphys.SpikeGroups = AnimalMetadata.ExtracellEphys.SpikeGroups;
    end
end

%Getting channel count - in a manner regardless of how user got spike
%groups... ie afterwards
channelcount = 0;
for gidx = 1:length(SessionMetadata.ExtracellEphys.SpikeGroups.groups)
    channelcount = channelcount+length(SessionMetadata.ExtracellEphys.SpikeGroups.groups{gidx});
end
SessionMetadata.ExtracellEphys.NumberOfChannels = channelcount;


%% Other modules like Virus, opto processing can happen here


%% Save all metadata
save(fullfile(basepath,[basename '.SessionMetadata.mat']),'SessionMetadata')


%% Handle and save sessionInfo (3/2019)
%should change this to use bz_editSessionInfo and use these as inputs to
%that
%maybe even abstract wrapper so field name and field value go into
%bz_editSession info as name-value pairs?

str = input('Write/Overwrite sessionInfo.mat file from SessionMetadata? (y/n): ','s');
if strcmp(str(1),'y')
    sessionInfo.session.name = basename;
    sessionInfo.session.path = basepath;
    sessionInfo.spikeGroups = SessionMetadata.ExtracellEphys.SpikeGroups;
    sessionInfo.nChannels = SessionMetadata.ExtracellEphys.NumberOfChannels;
    sessionInfo.channels = [0:SessionMetadata.ExtracellEphys.NumberOfChannels-1];
    sessionInfo.nBits = SessionMetadata.ExtracellEphys.Parameters.BitsPerSample;
    sessionInfo.rates.lfp = SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate;
    sessionInfo.rates.wideband = SessionMetadata.ExtracellEphys.Parameters.SampleRate;
    sessionInfo.rates.video = 0;%default
    sessionInfo.FileName = basename;
    sessionInfo.SampleTime = 1/SessionMetadata.ExtracellEphys.Parameters.SampleRate* 1e+6;%us apparently
    sessionInfo.nElecGps = [];
    for gidx = 1:(sessionInfo.spikeGroups.nGroups)%making elecgp... this is dumb
        tchans = sessionInfo.spikeGroups.groups{gidx};
        for cidx = 1:length(tchans)
            sessionInfo.ElecGp{gidx}.channel{cidx} = num2str(tchans(cidx));
        end
    end
    sessionInfo.HiPassFreq = [];
    sessionInfo.lfpSampleRate = SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate;
    sessionInfo.VoltageRange = SessionMetadata.ExtracellEphys.Parameters.VoltageRange;
    sessionInfo.Amplification = SessionMetadata.ExtracellEphys.Parameters.Amplification;
    sessionInfo.Offset = nan;
    d = AnimalMetadata.Surgery.Date;
    sessionInfo.Date = strcat(d(1:4),'-',d(5:6),'-',d(7:8));
    for gidx = 1:(sessionInfo.spikeGroups.nGroups)%making AnatGrps
        sessionInfo.AnatGrps(gidx).Channels = sessionInfo.spikeGroups.groups{gidx};
    end
    for gidx = 1:(sessionInfo.spikeGroups.nGroups)%making SpkGrps
        sessionInfo.SpkGrps(gidx).Channels = sessionInfo.spikeGroups.groups{gidx};
        sessionInfo.SpkGrps(gidx).nSamples = SessionMetadata.ExtracellEphys.Parameters.PointsPerWaveform;
        sessionInfo.SpkGrps(gidx).peakSample = SessionMetadata.ExtracellEphys.Parameters.PeakPointInWaveform;
        sessionInfo.SpkGrps(gidx).nFeatures = SessionMetadata.ExtracellEphys.Parameters.FeaturesPerWave;
    end
    %sessionInfo.Units = 
    sessionInfo.region = AnimalMetadata.ExtracellEphys.Channels.ChannelToAnatomyLookupTable.Table';
    sessionInfo.badchannels = SessionMetadata.ExtracellEphys.BadChannels;
%     sessionInfo.badshanks = SessionMetadata.ExtracellEphys.BadShanks;%not used now

    filename = fullfile(basepath,[basename,'.sessionInfo.mat']);
    save(filename,'sessionInfo'); 
    disp('sessionInfo successfully saved')
end

1;

function voltsperunit = VoltsPerUnit_Amplirec(basename,basepath)

%% Below for cases recorded in amplipex
metas = dir(fullfile(basepath,'*.meta'));
inis = dir(fullfile(basepath,'*.ini'));
if ~isempty(metas) | ~isempty(inis)
    bitdepth = 16;%just default
    bits = 2.^bitdepth;

    rangemax = [];
    rangemin = [];
    gain = [];

    if ~isempty(metas)%prefer to use .meta info, if not use .ini
        for a = 1:length(metas)
            fname = metas(a).name;
            if ~isempty(strfind(fname,basename))
                fname = fullfile(basepath,fname);
                fid=fopen(fname);
                tline= fgetl(fid);
                while ischar(tline)
                    try
                        if strcmp(tline(1:19),'Amplitude range max')
                            tline=tline(23:end);
                            rangemax(end+1)=str2num(tline);
                        end
                    end
                    try
                        if strcmp(tline(1:19),'Amplitude range min')
                            tline=tline(23:end);
                            rangemin(end+1)=str2num(tline);
                        end
                    end
                    try
                        if strcmp(tline(1:7),'Gain = ')
                            tline=tline(8:end);
                            gain(end+1) = str2num(tline);
                        end
                    end

                    tline= fgetl(fid);
                end
                fclose(fid);
            end
        end
    elseif ~isempty(inis)
        for a = 1:length(inis)
            fname = inis(a).name;
            if ~isempty(strfind(fname,basename))
                fname = fullfile(basepath,fname);
                fid=fopen(fname);
                tline= fgetl(fid);
                while ischar(tline)
                    try
                        if strcmp(tline(1:9),'rangeMax=')
                            tline=tline(10:end);
                            rangemax(end+1)=str2num(tline);
                        end
                    end
                    try
                        if strcmp(tline(1:9),'rangeMin=')
                            tline=tline(10:end);
                            rangemin(end+1)=str2num(tline);
                        end
                    end
                    try
                        if strcmp(tline(1:8),'auxGain=')
                            tline=tline(9:end);
                            gain(end+1) = str2num(tline);
                        end
                    end

                    tline= fgetl(fid);
                end
                fclose(fid);
            end
        end
    end

    if sum(abs(diff(rangemax))) | sum(abs(diff(rangemin))) 
        error('Ranges unequal between .ini files in this folder')
        return
    end
    if sum(abs(diff(gain)))
        error('Gains unequal between .ini files in this folder')
        return
    end

    if ~isempty(rangemax) & ~isempty(rangemin) %... if amplipex, otherwise do default which is intan
        totalrange = rangemax(1)-rangemin(1);
        voltsperunit = totalrange / bits / gain(1);
    end
end    
