function SessionMetadata = bz_SetSessionMetadata(basepath,basename)

if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
if ~exist('basename','var')
    [~,basename] = fileparts(basepath);
elseif isempty(basename)
    [~,basename] = fileparts(basepath);
end

%% Manual entries
BadShanks = [];
BadChannels = [];
NumberOfTurnsSinceSurgery = [];
ChannelNotes = {''};

%% First step, look up one directory to look for AnimalMetadata and copy it 
% into this folder - this will give crucial info for subsequent work here
% Make an initial version of output from this, overwrite with actual
% metadata from recording system
% 
supradir = fileparts(basepath);
d = dir(fullfile(supradir,'*_AnimalMetadata.mat'));
upperamdpath = fullfile(supradir,d(1).name);
thisamdpath = fullfile(basepath,[basename,'_AnimalMetadata.mat']);
copyfile(upperamdpath,thisamdpath);
load(thisamdpath);%will now have the variable AnimalMetadata

%copy xml file from animal folder to here
d = dir(fullfile(supradir,'*.xml'));
upperxmlpath = fullfile(supradir,d(1).name);
thisxmlpath = fullfile(basepath,[basename,'.xml']);
copyfile(upperxmlpath,thisxmlpath);

% set basic info based on AnimalMetadata
SessionMetadata.Basepath = basepath;
SessionMetadata.Basename = basename;
SessionMetadata.AnimalBasepath = AnimalMetadata.AnimalBasepath;
SessionMetadata.Weight = [];
SessionMetadata.RecordingParameters = AnimalMetadata.RecordingParameterDefaults;

SessionMetadata.Anatomy.TargetRegions = AnimalMetadata.Probes.TargetRegions;
SessionMetadata.Anatomy.NumberOfTurnsSinceSurgery = [0 0];
SessionMetadata.Anatomy.ProbeDepths = SessionMetadata.Anatomy.NumberOfTurnsSinceSurgery .* AnimalMetadata.Probes.UmPerScrewTurn;
SessionMetadata.Anatomy.ProbeCenterCoordinates = [];
SessionMetadata.Anatomy.ChannelCenterCoordinates = [];

SessionMetadata.Probes.ProbeLayoutFilenames = AnimalMetadata.Probes.ProbeLayoutFilenames;
SessionMetadata.Probes.NumChansPerProbe = AnimalMetadata.Probes.NumChansPerProbe;
SessionMetadata.Probes.ProbeSpikeGroupLayoutSuperficialToDeep = AnimalMetadata.Probes.ProbeSpikeGroupLayoutSuperficialToDeep;
SessionMetadata.Probes.BadShanks = [];%manual... can be used to make badchannels

SessionMetadata.Channels.ChannelToProbeLookupTable_AsPlugged = AnimalMetadata.Channels.ChannelToProbeLookupTable_AsPlugged;
SessionMetadata.Channels.ChannelToGroupLookupTable_AsPlugged = AnimalMetadata.Channels.ChannelToGroupLookupTable_AsPlugged;
SessionMetadata.Channels.ChannelToAnatomyLookupTable_AsPlugged = AnimalMetadata.Channels.ChannelToAnatomyLookupTable_AsPlugged;
SessionMetadata.Channels.BadChannels = BadChannels;%manual list in addition to channels that can be derived from .Probes.BadShanks
SessionMetadata.Channels.ChannelNotes = ChannelNotes;%manual list in addition to channels that can be derived from .Probes.BadShanks

%what else to add from AnimalMetadata... not sure if would rather be
%redundant or not

%% Read recording system-based metadata from either .rhd or .meta (Intan or Amplirec) 
% Check for compatibility with AnimalMetadata, overwrite using this and
% warn the user of the conflict.

d = dir(fullfile(basepath,'*.meta'));
SessionMetadata.RecordingFiles.Names = {};
SessionMetadata.RecordingFiles.Bytes = [];
SessionMetadata.RecordingFiles.Seconds = [];

if ~isempty(d);%if from an amplipex
    SessionMetadata.RecordingSystem = 'Amplipex';
    for idx = 1:length(d);
        SessionMetadata.RecordingFiles.Names{idx} = strcat(d(idx).name(1:end-4),'dat');
        t = dir(fullfile(basepath,SessionMetadata.RecordingFiles.Names{idx}));
        if ~isempty(t) %if original .dat still around
            SessionMetadata.RecordingFiles.Bytes(idx) = t.bytes;
        else%if no .dat around
            disp([SessionMetadata.RecordingFiles.Names{idx} ' not found']);
            SessionMetadata.RecordingFiles.Bytes(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'filebytes'));
        end
        
        NChannels(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'NumChans'));
        SamplingRate(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'SamplingRate'));
        Amplification(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'Gain'));
    end
    if ~isempty(NChannels) %if some files were found
        %Channels
        if sum(abs(diff(NChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers of channels - per .meta files.  Will use count from first recording.')
        else%if all the same
            SessionMetadata.RecordingParameters.NumberOfChannels = NChannels(1);
            if NChannels(1) ~= SessionMetadata.RecordingParameters.NumberOfChannels
                warning('The number of channels recorded is different from that in AnimalMetadata.  Original XML File invalid');
            end
        end

        %Sampling Rate
        if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different SamplingRates - per .meta files.  Will use number from first recording.')
        else%if all the same
            SessionMetadata.RecordingParameters.SampleRate = SamplingRate(1);
            if SamplingRate(1) ~= SessionMetadata.RecordingParameters.SampleRate
                warning('The Sampling Rates recorded is different from that in AnimalMetadata.  Original XML File invalid');
            end
        end
        
        %Amplification
        if sum(abs(diff(Amplification)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different Amplifications - per .meta files.  Will use number from first recording.')
        else%if all the same
            SessionMetadata.RecordingParameters.Amplification = Amplification(1);
            if Amplification(1) ~= SessionMetadata.RecordingParameters.Amplification
                warning('The number of channels recorded is different from that in AnimalMetadata.  Original XML File invalid');
            end
        end
        
        %Not changing other RecordingParameters

        % getting bytes and seconds per file
        SessionMetadata.VoltsPerUnit = VoltsPerUnit_Amplirec(basename,basepath);%calculate from .metas/.inis
    end
else %basically if intan
    NAmpChannels = [];
    NAuxChannels = [];
    NDigitalChannels = [];
    SamplingRate = [];
    d = dir(basepath);
    for didx = 1:length(d)
       if d(didx).isdir
           t = dir(fullfile(basepath,d(didx).name,'info.rhd'));
           if ~isempty(t)
               SessionMetadata.RecordingSystem = 'Intan';
               rhd = fullfile(basepath,d(didx).name,'info.rhd');
%                ampdats{end+1} = fullfile(basepath,d(didx).name,'amplifier.dat');
%                recordingdiridxs = didx;
               SessionMetadata.RecordingFiles.Names{end+1} = [d(didx).name '.dat'];

               [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
               board_dig_in_channels, supply_voltage_channels, frequency_parameters ] =...
               read_Intan_RHD2000_file(fullfile(basepath,d(didx).name),'info.rhd');
           
               NAmpChannels(end+1) = length(amplifier_channels);
               NAuxChannels(end+1) = length(aux_input_channels);
               NDigitalChannels(end+1) = length(board_dig_in_channels);
               SamplingRate(end+1) = frequency_parameters.amplifier_sample_rate;
               
               t = dir(fullfile(basepath,d(didx).name,'amplifier.dat'));
               if ~isempty(t) %if original .dat still around
                   SessionMetadata.RecordingFiles.Bytes(end+1) = t.bytes;
               else%if no .dat around
                   disp([SessionMetadata.RecordingFiles.Names{end} ' not found']);
                   timedat = dir(fullfile(basepath,d(didx).name,'time.dat'));
                   SessionMetadata.RecordingFiles.Bytes(end+1) = timedat(1).bytes/2*NAmpChannels(end);
               end

           end
        end
    end
    if ~isempty(NAmpChannels) %if some files were found
        %channels
        if sum(abs(diff(NAmpChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers of channels - per info.rhd files.  Will use count from first recording.')
        else%if all the same
            SessionMetadata.RecordingParameters.NumberOfChannels = NAmpChannels(1);
            if NAmpChannels(1) ~= SessionMetadata.RecordingParameters.NumberOfChannels
                warning('The number of channels recorded is different from that in AnimalMetadata.  Original XML File invalid');
            end
        end

        SessionData.Channels.NumAuxiliaryChannels = NAuxChannels(1);
        if sum(abs(diff(NAuxChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers auxiliary of channels - per info.rhd files.  Will use count from first recording.')
        end

        SessionData.Channels.NumDigitalChannels = NAuxChannels(1);
        if sum(abs(diff(NDigitalChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers digital of channels - per info.rhd files.  Will use count from first recording.')
        end

        %sampling rate
        if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different SamplingRates - per info.rhd files.  Will use number from first recording.')
        else%if all the same
            SessionMetadata.RecordingParameters.SampleRate = SamplingRate(1);
            if SamplingRate(1) ~= SessionMetadata.RecordingParameters.SampleRate
                warning('The number of Sampling Rates recorded is different from that in AnimalMetadata.  Original XML File invalid');
            end
        end
        %Not changing other RecordingParameters

        % getting bytes and seconds per file
        SessionMetadata.VoltsPerUnit = 0.0000002;%default for Intan    
    end
end

SessionMetadata.RecordingFiles.Seconds = SessionMetadata.RecordingFiles.Bytes/SessionMetadata.RecordingParameters.NumberOfChannels/2/SessionMetadata.RecordingParameters.SampleRate;

save(fullfile(basepath,[basename '_SessionMetadata.mat']),'SessionMetadata')

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
