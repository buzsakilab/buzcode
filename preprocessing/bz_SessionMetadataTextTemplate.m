function SessionMetadata = bz_SessionMetadataTextTemplate(basepath)
% This .m file is a generic template and will be copied to animal
% folders/session folders.  Edits here will change what is automatically
% copied to animal folders when bz_EditAnimalMetadata.m is called.
% 
% Developers should edit bz_SessionMetadataTextTemplate.m to change how
% SessionMetadata is created.
% 
% Brendon Watson 2017


if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

%% HUMAN INPUT BELOW
% Surgery and Animal metadata
SessionMetadata.Animal.WeightGrams = [250];

% Extracell Ephys metadata
SessionMetadata.ExtracellEphys.NumberOfTurnsSinceSurgery = [0];%vector, one entry per probe
SessionMetadata.ExtracellEphys.Probes.PluggingOrder = [];%vector, one entry per probe. blank defaults to animal plugging order
SessionMetadata.ExtracellEphys.BadShanks = [];% vector for this recording. base 1
     % These bad shanks will be used to populate bad channels
SessionMetadata.ExtracellEphys.BadChannels = [];% vector for this recording. base 0
SessionMetadata.ExtracellEphys.ChannelNotes = {''};

SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate = 1250;%assumed default
SessionMetadata.ExtracellEphys.Parameters.PointsPerWaveform = 32;%default
SessionMetadata.ExtracellEphys.Parameters.PeakPointInWaveform = 16;%default
SessionMetadata.ExtracellEphys.Parameters.FeaturesPerWave = 4;%default

% Optogenetics metadata
FiberNum = 1;%copy more of these blocks - one per probe
SessionMetadata.Optogenetics.Fibers(FiberNum).StimulusRecordingChannelDatType = 'AnalogIn';%'AnalogIn' for Intan, 'Main' if fused with other dat
SessionMetadata.Optogenetics.Fibers(FiberNum).StimulusRecordingChannelNumber = [];%
%?else

% Behavior metadata
BehaviorSessionNumber = 1;
SessionMetadata.BehaviorEpisodes.(BehaviorSessionNumber).StartStopSeconds = [];%startstop pair, in seconds
SessionMetadata.BehaviorEpisodes(BehaviorSessionNumber).DatFileIndex = [];%which dat files are included in this Episode
SessionMetadata.BehaviorEpisodes(BehaviorSessionNumber).BehaviorType = '';%verbal description of behavior type, ie "LinearTrack"
SessionMetadata.BehaviorEpisodes(BehaviorSessionNumber).BehaviorNotes = '';
SessionMetadata.BehaviorEpisodes.(BehaviorSessionNumber).Motion.MotionTrackerType = '';%Optitrack or Basler etc
SessionMetadata.BehaviorEpisodes.(BehaviorSessionNumber).Motion.CameraFrameChannelType = 'DigitalIn';%'DigitalIn' for Intan, 'Main' if fused with other dat
SessionMetadata.BehaviorEpisodes.(BehaviorSessionNumber).Motion.CameraFrameChannelNumber = [];%channel number where frame times are recorded
SessionMetadata.BehaviorEpisodes(BehaviorSessionNumber).Location = '';

SleepSessionNumber = 1;
SessionMetadata.SleepEpisodes.(SleepSessionNumber).StartStopSeconds = [];%startstop pair, in seconds
SessionMetadata.SleepEpisodes(SleepSessionNumber).DatFileIndex = [];%which dat files are included in this Episode
SessionMetadata.SleepEpisodes(SleepSessionNumber).BehaviorNotes = '';
SessionMetadata.SleepEpisodes.(SleepSessionNumber).Motion.MotionTrackerType = '';%Optitrack or Basler etc
SessionMetadata.SleepEpisodes.(SleepSessionNumber).Motion.CameraFrameChannelType = 'DigitalIn';%'DigitalIn' for Intan, 'Main' if fused with other dat
SessionMetadata.SleepEpisodes.(SleepSessionNumber).Motion.CameraFrameChannelNumber = [];%channel number where frame times are recorded
SessionMetadata.SleepEpisodes(SleepSessionNumber).Location = '';

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
d1 = dir(fullfile(basepath,'*.AnimalMetadata.mat'));
d2 = dir(fullfile(supradir,'*.AnimalMetadata.mat'));
guessamdpath = [];
finalamdpath = fullfile(basepath,[basename,'.AnimalMetadata.mat']);
if ~isempty(d1)%if already one in this path
    guessamdpath = fullfile(basepath,d1(1).name);
elseif ~isempty(d2)%if one in directory above
    guessamdpath = fullfile(supradir,d2(1).name);
else%if cannot be found, ask user 
    disp('No local AnimalMetadata.mat found');
end
[FileName,PathName] = uigetfile('.AnimalMetadata.mat','Find AnimalMetaData.mat',guessamdpath);
inputamdpath = fullfile(PathName,FileName);

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
    SessionMetadata.ExtracellEphys.Files.Names = {};
    SessionMetadata.ExtracellEphys.Files.Bytes = [];
    SessionMetadata.ExtracellEphys.Files.Seconds = [];

    %first check if amplipex (.metas present) or intan (rhd.infos
    %present)... or neither
    recsys = '';
    
    dmeta = dir(fullfile(basepath,'*.meta'));
    if ~isempty(dmeta)
        recsys = 'Amplipex';
    end
    
    drhd = dir(basepath);
    for didx = 1:length(drhd)
       if drhd(didx).isdir
           t = dir(fullfile(basepath,drhd(didx).name,'info.rhd'));
           if ~isempty(t)%only use folders with info.rhd inside
                recsys = 'Intan';
           end
       end
    end
    
    
    switch recsys
        case 'Amplipex'%if from an amplipex
            d = dmeta;
            SessionMetadata.ExtracellEphys.RecordingSystem = 'Amplipex';
            for idx = 1:length(d);
                SessionMetadata.ExtracellEphys.Files.Names{idx} = d(idx).name(1:end-5);
                tdat = strcat(d(idx).name(1:end-5),'.dat');
                t = dir(fullfile(basepath,tdat));
                if ~isempty(t) %if original .dat still around
                    SessionMetadata.ExtracellEphys.Files.Bytes(idx) = t.bytes;
                else%if no .dat around (ie anymore), get bytes from .meta 
                    disp([tdat ' not found']);
                    SessionMetadata.ExtracellEphys.Files.Bytes{idx} = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'filebytes'));
                end
                NChannels(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'NumChans'));
                SamplingRate(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'SamplingRate'));
                Amplification(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'Gain'));
            end
            if ~isempty(NChannels) %if some files were found
                %Channels
                if sum(abs(diff(NChannels)))%if different numbers of channels differed across recordings in session, error out
                    warning('Rercordings contain different numbers of channels - per .meta files.  Will use count from first recording.')
                end
                SessionMetadata.ExtracellEphys.Parameters.NumberOfChannels = NChannels(1);

                %Sampling Rate
                if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
                    warning('Rercordings contain different SamplingRates - per .meta files.  Will use number from first recording.')
                end
                SessionMetadata.ExtracellEphys.Parameters.SampleRate = SamplingRate(1);

                %Amplification
                if sum(abs(diff(Amplification)))%if different numbers of channels differed across recordings in session, error out
                    warning('Rercordings contain different Amplifications - per .meta files.  Will use number from first recording.')
                end
                SessionMetadata.ExtracellEphys.Parameters.Amplification = Amplification(1);

                SessionMetadata.ExtracellEphys.Parameters.VoltsPerUnit = VoltsPerUnit_Amplirec(basename,basepath);%calculate from .metas/.inis
                SessionMetadata.ExtracellEphys.Parameters.BitsPerSample = 16;%amplirec default
                SessionMetadata.ExtracellEphys.Parameters.VoltageRange = 10;%not used except to make xml
            end
        case 'Intan' %if not amplipex, look for intan
            d = drhd;
            NAmpChannels = [];
            NAuxChannels = [];
            NDigitalChannels = [];
            SamplingRate = [];

            for didx = 1:length(d)
                if d(didx).isdir
                    t = dir(fullfile(basepath,d(didx).name,'info.rhd'));
                    if ~isempty(t)%only use folders with info.rhd inside

                       SessionMetadata.ExtracellEphys.RecordingSystem = 'Intan';
                       rhd = fullfile(basepath,d(didx).name,'info.rhd');
        %                ampdats{end+1} = fullfile(basepath,d(didx).name,'amplifier.dat');
        %                recordingdiridxs = didx;
                       SessionMetadata.ExtracellEphys.Files.Names{end+1} = d(didx).name;

                       % read the info.rhd file for each individual recorded
                       % file/folder
                       try
                          [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
                              board_dig_in_channels, supply_voltage_channels, frequency_parameters ] =...
                              read_Intan_RHD2000_file(fullfile(basepath,d(didx).name),'info.rhd');%function from Intan
                          NAmpChannels(end+1) = length(amplifier_channels);
                          NAuxChannels(end+1) = length(aux_input_channels);
                          NDigitalChannels(end+1) = length(board_dig_in_channels);
                          SamplingRate(end+1) = frequency_parameters.amplifier_sample_rate;

                       catch%if read error, set values based on prior
                           if ~isempty(NAmpChannels)
    %                            if ~isnan(NAmpChannels(1));
                                  NAmpChannels(end+1) = NAmpChannels(1);
    %                            else
    %                               NAmpChannels(end+1) = SessionMetadata.AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal;
    %                            end
                           else
                               NAmpChannels(end+1) = nan;
                           end
                           if ~isempty(NAuxChannels)
                               NAuxChannels(end+1) = NAuxChannels(1);
                           else
                               NAuxChannels(end+1) = nan;
                           end
                           if ~isempty(NDigitalChannels)
                               NDigitalChannels(end+1) = NDigitalChannels(1);
                           else
                               NDigitalChannels(end+1) = nan;
                           end
                           if ~isempty(SamplingRate)
                               SamplingRate(end+1) = SamplingRate(1);
                           else
                               SamplingRate(end+1) = nan;
                           end
                       end

                       t = dir(fullfile(basepath,d(didx).name,'amplifier.dat'));
                       if ~isempty(t) %if original .dat still around
                           SessionMetadata.ExtracellEphys.Files.Bytes(end+1) = t.bytes;
                       else%if no .dat around
                           disp([SessionMetadata.ExtracellEphys.Files.Names{end} ' not found']);
                           timedat = dir(fullfile(basepath,d(didx).name,'time.dat'));
                           if ~isnan(NAmpChannels(end));
                              SessionMetadata.ExtracellEphys.Files.Bytes(end+1) = timedat(1).bytes/2*NAmpChannels(end);
                           else
                              SessionMetadata.ExtracellEphys.Files.Bytes(end+1) = timedat(1).bytes/2*SessionMetadata.AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal;
                           end
                       end

                   end
                end
            end
            %handle errors created by improper read of rhd.info files
            if isnan(sum(abs(diff(NAmpChannels))))
                warning('At least one read error on rhd.info files.  Using values from other files to populate.')
                ok = find(~isnan(NAmpChannels),1,'first');
                if isempty(ok)
                    warning('No good reads on NumberofAmpChannels. Assuming based on probe maps')
                    NAmpChannels(isnan(NAmpChannels)) = SessionMetadata.AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal
                else
                    NAmpChannels(isnan(NAmpChannels)) = NAmpChannels(ok);
                end
            end
            if isnan(sum(abs(diff(NAuxChannels))))
                warning('At least one read error on rhd.info files.  Using values from other files to populate.')
                ok = find(~isnan(NAuxChannels),1,'first');
                if isempty(ok)
                    error('No good reads on NAuxChannels, assuming 0')
                    NAuxChannels(isnan(NAuxChannels)) = 0;
                else
                    NAuxChannels(isnan(NAuxChannels)) = NAuxChannels(ok);
                end
            end
            if isnan(sum(abs(diff(NDigitalChannels))))
                warning('At least one read error on rhd.info files.  Using values from other files to populate.')
                ok = find(~isnan(NDigitalChannels),1,'first');
                if isempty(ok)
                    error('No good reads on NDigitalChannels, assuming 0')
                    NDigitalChannels(isnan(NDigitalChannels)) = 0;
                else
                    NDigitalChannels(isnan(NDigitalChannels)) = NDigitalChannels(ok);
                end
            end
            if isnan(sum(abs(diff(SamplingRate))))
                warning('At least one read error on rhd.info files.  Using values from other files to populate.')
                ok = find(~isnan(SamplingRate),1,'first');
                if isempty(ok)
                    error('No good reads on SamplingRate, taking from AnimalMetadata')
                    SamplingRate(isnan(SamplingRate)) = SessionMetadata.AnimalMetadata.ExtracellEphys.Parameters.SampleRate;
                else
                    SamplingRate(isnan(SamplingRate)) = SamplingRate(ok);
                end
            end

            %handle errors where different parameters read for different files (rather than uniform)            
            if sum(abs(diff(NAmpChannels)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different numbers of channels - per info.rhd files.  Will use count from first recording.')
            end
            SessionMetadata.ExtracellEphys.Parameters.NumberOfChannels = NAmpChannels(1);

            if sum(abs(diff(NAuxChannels)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different numbers auxiliary of channels - per info.rhd files.  Will use count from first recording.')
            end
            SessionMetadata.ExtracellEphys.NumAuxiliaryChannels = NAuxChannels(1);

            if sum(abs(diff(NDigitalChannels)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different numbers digital of channels - per info.rhd files.  Will use count from first recording.')
            end
            SessionMetadata.ExtracellEphys.NumDigitalChannels = NAuxChannels(1);

            %sampling rate
            if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different SamplingRates - per info.rhd files.  Will use number from first recording.')
            end
            SessionMetadata.ExtracellEphys.Parameters.SampleRate = SamplingRate(1);

            % Not gather-able from intan info.rhd files
            SessionMetadata.ExtracellEphys.Parameters.Amplification = 1;%digitized on chip, let's say 1
            SessionMetadata.ExtracellEphys.Parameters.VoltsPerUnit = 0.0000002;%intan default                
            SessionMetadata.ExtracellEphys.Parameters.BitsPerSample = 16;%intan default
            SessionMetadata.ExtracellEphys.Parameters.VoltageRange = 10;%not used except to make xml

            %Now that everything is in order, calculate seconds from bytes
            SessionMetadata.ExtracellEphys.Files.Seconds = SessionMetadata.ExtracellEphys.Files.Bytes/SessionMetadata.ExtracellEphys.Parameters.NumberOfChannels/2/SessionMetadata.ExtracellEphys.Parameters.SampleRate;        
        case '' %if no known recording system used
            warning('No recording system found, using defaults.  May manually set values in guts of this .m file')
            SessionMetadata.ExtracellEphys.RecordingSystem = 'Unknown';
            SessionMetadata.ExtracellEphys.Files.Names = {};
            SessionMetadata.ExtracellEphys.Files.Bytes = [];
            SessionMetadata.ExtracellEphys.Files.Seconds = [];%can manually write in here if desired

            SessionMetadata.ExtracellEphys.Parameters.NumberOfChannels = SessionMetadata.AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal;

            SessionMetadata.ExtracellEphys.Parameters.SampleRate = 20000;%Lab default
            SessionMetadata.ExtracellEphys.Parameters.Amplification = 1;%Intan digitized on chip, let's say 1
            SessionMetadata.ExtracellEphys.Parameters.VoltsPerUnit = 0.0000002;%Intan default                
            SessionMetadata.ExtracellEphys.Parameters.BitsPerSample = 16;%Intan default
            SessionMetadata.ExtracellEphys.Parameters.VoltageRange = 10;%not used except to make xml
    end

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
    if isempty(SessionMetadata.ExtracellEphys.BadChannels)
        if ~isempty(SessionMetadata.ExtracellEphys.BadShanks)
            chanpergp = SessionMetadata.AnimalMetadata.ExtracellEphys.Probes.ProbeSpikeGroupLayoutSuperficialToDeep;
            badchans = [];
            for sidx = 1:length(SessionMetadata.ExtracellEphys.BadShanks)
                tshank = SessionMetadata.ExtracellEphys.BadShanks;
                badchans = cat(1,badchans,chanpergp{tshank});
            end
            SessionMetadata.ExtracellEphys.BadChannels = badchans;
        end
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
