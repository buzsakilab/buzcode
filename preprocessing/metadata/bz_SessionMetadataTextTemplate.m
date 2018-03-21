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


%% HUMAN INPUT BELOW
% Preprocessing specifiers
SessionMetadata.Preprocess.AnimalMetadataSource = {'AnimalFolderAbove'};

% Surgery and Animal metadata
SessionMetadata.Animal.WeightGrams = [];

% Extracell Ephys metadata
SessionMetadata.ExtracellEphys.NumberOfTurnsSinceSurgery = [0];%vector, one entry per probe
SessionMetadata.ExtracellEphys.Probes.PluggingOrder = [];%vector, one entry per probe. blank defaults to animal plugging order
SessionMetadata.ExtracellEphys.BadShanks = [];% vector for this recording. base 1
     % These bad shanks will be used to populate bad channels
SessionMetadata.ExtracellEphys.BadChannels = [];% vector for this recording. base 0
SessionMetadata.ExtracellEphys.ChannelNotes = {''};
SessionMetadata.ExtracellEphys.SpikeGroups = {'FromAnimalMetaData','FromXML'}; %Pick one, delete the other.

SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate = 1250;%assumed default
SessionMetadata.ExtracellEphys.Parameters.PointsPerWaveform = 32;%default
SessionMetadata.ExtracellEphys.Parameters.PeakPointInWaveform = 16;%default
SessionMetadata.ExtracellEphys.Parameters.FeaturesPerWave = 4;%default

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
if strcmp(SessionMetadata.Preprocess.AnimalMetadataSource{1},'AnimalFolderAbove')
    d2 = dir(fullfile(supradir,'*.AnimalMetadata.mat'));
    inputamdpath = fullfile(supradir,d2(1).name);
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

    if isempty(SessionMetadata.ExtracellEphys.Probes.PluggingOrder)
        SessionMetadata.ExtracellEphys.Probes.PluggingOrder = SessionMetadata.AnimalMetadata.ExtracellEphys.Probes.PluggingOrder;
    end
    
    % Read recording system-based metadata from either .rhd or .meta (Intan or Amplirec)
    % Check for compatibility with AnimalMetadata, overwrite using this and
    % warn the user of the conflict.
    SessionMetadata.ExtracelEphys = DatInfoMake(basepath);

    %now things regardless of input system

    %Check for conflicts between original animal-based ephys params and
    %those based on files in this particular recording.
    aparams = SessionMetadata.AnimalMetadata.ExtracellEphys.Parameters;
    sparams = SessionMetadata.ExtracellEphys.Parameters;
    fn = fieldnames(aparams);
    divergentfields = {};
    for fidx = 1:length(fn);
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
    
    %Get badchannels from badshanks... ie all channels on bad shanks are
    %bad
    % Add bad channels to channels derived from bad channels
%     if isempty(SessionMetadata.ExtracellEphys.BadChannels)
        if ~isempty(SessionMetadata.ExtracellEphys.BadShanks)
            chanpergp = SessionMetadata.AnimalMetadata.ExtracellEphys.Probes.ProbeSpikeGroupLayoutSuperficialToDeep;
            badchans = [];
            for sidx = 1:length(SessionMetadata.ExtracellEphys.BadShanks)
                tshank = SessionMetadata.ExtracellEphys.BadShanks;
                badchans = cat(1,badchans,chanpergp{tshank});
            end
            SessionMetadata.ExtracellEphys.BadChannels = cat(1,SessionMetadata.ExtracellEphys.BadChannels(:),badchans(:));
        end
%     end 
    
    %Get Spikegroups
    if strcmp(SessionMetadata.ExtracellEphys.SpikeGroups,'FromXML'); 
        display('Getting spikegroups from the .xml...')
        %Load the XML for spike groups
        xmlname = fullfile(basepath,[basename,'.xml']);
        if ~exist(xmlname,'file')
            [FileName,PathName] = uigetfile('.xml','Find the .xml to load SpikeGroups from',basepath);
            xmlname = fullfile(PathName,FileName);
        end
        xmlparms = LoadParameters(xmlname);
        SessionMetadata.ExtracellEphys.SpikeGroups = xmlparms.SpkGrps;
    end
end

%% Other modules like Virus, opto processing can happen here


%% Save all metadata
save(fullfile(basepath,[basename '.SessionMetadata.mat']),'SessionMetadata')

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
