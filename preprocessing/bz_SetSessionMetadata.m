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
    
%% First step, for AnimalMetadata to open it.  Copy it to basepath if it's not already there. 
% into this folder - this will give crucial info for subsequent work here
% Make an initial version of output from this, overwrite with actual
% metadata from recording system
% 
supradir = fileparts(basepath);
d1 = dir(fullfile(basepath,'*_AnimalMetadata.mat'));
d2 = dir(fullfile(supradir,'*_AnimalMetadata.mat'));
if ~isempty(d1)%if already one in this path
    thisamdpath = fullfile(basepath,d1(1).name);
elseif ~isempty(d2)%if one in directory above
    upperamdpath = fullfile(supradir,d2(1).name);
    thisamdpath = fullfile(basepath,[basename,'_AnimalMetadata.mat']);
    copyfile(upperamdpath,thisamdpath);
else%if cannot be found, ask user 
    disp('Find the AnimalMetadata.mat file for this animal');
    [FileName,PathName] = uigetfile('_AnimalMetadata.mat','Find AnimalMetaData.mat');
    thatamdpath = fullfile(FileName,PathName);
    copyfile(thatamdpath,thisamdpath);
end
load(thisamdpath);
clear d1 d2 supradir upperamdpath FileName PathName thatamdpath thisamdpath

%% Value setting - humans must do this.  Will be asked to edit a _NoteText.m file
notesname = [basename,'_SessionNotesText'];
notesfullpath = fullfile(basepath,[basename,'_SessionNotesText.m']);
if ~exist(fullfile(basepath,notesfullpath),'file')
    w = which('bz_SessionNotesTemplate.m');% copy an example header here to edit
    copyfile(w,notesfullpath);
end

edit(notesfullpath)
prompt = 'Push any key in this window when done editing the NoteText file ';
str = input(prompt,'s');
eval([notesname '(basepath,basename,AnimalMetadata);']);%save _AnimalNotes.mat to disk

load(fullfile(basepath,[basename '_SessionNotes.mat']))%load AnimalNotes
SessionMetadata = SessionNotes;
SessionMetadata.AnimalMetadata = AnimalMetadata;%will now have the variable AnimalMetadata

clear SessionNotes AnimalMetadata

%%
if SessionMetadata.AnimalMetadata.Modules.ExtracellEphys
    SessionMetadata.Probes.NumberOfTurnsSinceSurgery = [0 0];
    SessionMetadata.Probes.ProbeDepths = SessionMetadata.ExtraEphys.NumberOfTurnsSinceSurgery .* SessionMetadata.AnimalMetadata.Probes.UmPerScrewTurn;
    SessionMetadata.Probes.ProbeCenterCoordinates = [];%someone should do this
    SessionMetadata.Probes.ChannelCenterCoordinates = [];%someone should do this

    % Read recording system-based metadata from either .rhd or .meta (Intan or Amplirec) 
    % Check for compatibility with AnimalMetadata, overwrite using this and
    % warn the user of the conflict.

    SessionMetadata.ExtraEphys.Files.Names = {};
    SessionMetadata.ExtraEphys.Files.Bytes = [];
    SessionMetadata.ExtraEphys.Files.Seconds = [];

    d = dir(fullfile(basepath,'*.meta'));
    if ~isempty(d);%if from an amplipex
        SessionMetadata.ExtraEphys.RecordingSystem = 'Amplipex';
        for idx = 1:length(d);
            SessionMetadata.ExtraEphys.Files.Names{idx} = d(idx).name(1:end-5);
            tdat = strcat(d(idx).name(1:end-5),'.dat');
            t = dir(fullfile(basepath,tdat));
            if ~isempty(t) %if original .dat still around
                SessionMetadata.ExtraEphys.Files.Bytes(idx) = t.bytes;
            else%if no .dat around (ie anymore), get bytes from .meta 
                disp([tdat ' not found']);
                SessionMetadata.ExtraEphys.Files.Bytes{idx} = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'filebytes'));
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
            SessionMetadata.ExtraEphys.Parameters.NumberOfChannels = NChannels(1);
            
            %Sampling Rate
            if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different SamplingRates - per .meta files.  Will use number from first recording.')
            end
            SessionMetadata.ExtraEphys.Parameters.SampleRate = SamplingRate(1);
            
            %Amplification
            if sum(abs(diff(Amplification)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different Amplifications - per .meta files.  Will use number from first recording.')
            end
            SessionMetadata.ExtraEphys.Parameters.Amplification = Amplification(1);
            
            SessionMetadata.ExtraEphys.Parameters.VoltsPerUnit = VoltsPerUnit_Amplirec(basename,basepath);%calculate from .metas/.inis
            SessionMetadata.ExtraEphys.Parameters.BitsPerSample = 16;%amplirec default
            SessionMetadata.ExtraEphys.Parameters.VoltageRange = 10;%not used except to make xml
        end
    else %if not amplipex, look for intan
        NAmpChannels = [];
        NAuxChannels = [];
        NDigitalChannels = [];
        SamplingRate = [];
        d = dir(basepath);
        for didx = 1:length(d)
           if d(didx).isdir
               t = dir(fullfile(basepath,d(didx).name,'info.rhd'));
               if ~isempty(t)%only use folders with info.rhd inside
                   SessionMetadata.ExtraEphys.RecordingSystem = 'Intan';
                   rhd = fullfile(basepath,d(didx).name,'info.rhd');
    %                ampdats{end+1} = fullfile(basepath,d(didx).name,'amplifier.dat');
    %                recordingdiridxs = didx;
                   SessionMetadata.ExtraEphys.Files.Names{end+1} = d(didx).name;

                   % read the info.rhd file for each individual recorded
                   % file/folder
                   [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
                   board_dig_in_channels, supply_voltage_channels, frequency_parameters ] =...
                   read_Intan_RHD2000_file(fullfile(basepath,d(didx).name),'info.rhd');%function from Intan

                   NAmpChannels(end+1) = length(amplifier_channels);
                   NAuxChannels(end+1) = length(aux_input_channels);
                   NDigitalChannels(end+1) = length(board_dig_in_channels);
                   SamplingRate(end+1) = frequency_parameters.amplifier_sample_rate;

                   t = dir(fullfile(basepath,d(didx).name,'amplifier.dat'));
                   if ~isempty(t) %if original .dat still around
                       SessionMetadata.ExtraEphys.Files.Bytes(end+1) = t.bytes;
                   else%if no .dat around
                       disp([SessionMetadata.ExtraEphys.Files.Names{end} ' not found']);
                       timedat = dir(fullfile(basepath,d(didx).name,'time.dat'));
                       SessionMetadata.ExtraEphys.Files.Bytes(end+1) = timedat(1).bytes/2*NAmpChannels(end);
                   end

               end
            end
        end
        if ~isempty(NAmpChannels) %if some files were found
            %channels
            if sum(abs(diff(NAmpChannels)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different numbers of channels - per info.rhd files.  Will use count from first recording.')
            end
            SessionMetadata.ExtraEphys.Parameters.NumberOfChannels = NAmpChannels(1);

            if sum(abs(diff(NAuxChannels)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different numbers auxiliary of channels - per info.rhd files.  Will use count from first recording.')
            end
            SessionData.ExtraEphys.NumAuxiliaryChannels = NAuxChannels(1);
            
            if sum(abs(diff(NDigitalChannels)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different numbers digital of channels - per info.rhd files.  Will use count from first recording.')
            end
            SessionData.ExtraEphys.NumDigitalChannels = NAuxChannels(1);

            %sampling rate
            if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different SamplingRates - per info.rhd files.  Will use number from first recording.')
            end
            SessionMetadata.ExtraEphys.Parameters.SampleRate = SamplingRate(1);

            % getting bytes and seconds per file
            SessionMetadata.ExtraEphys.Parameters.Amplification = 1;%digitized on chip, let's say 1
            SessionMetadata.ExtraEphys.Parameters.VoltsPerUnit = 0.0000002;%intan default                
            SessionMetadata.ExtraEphys.Parameters.BitsPerSample = 16;%intan default
            SessionMetadata.ExtraEphys.Parameters.VoltageRange = 10;%not used except to make xml
        end
    end

    SessionMetadata.ExtraEphys.Files.Seconds = SessionMetadata.ExtraEphys.Files.Bytes/SessionMetadata.ExtraEphys.Parameters.NumberOfChannels/2/SessionMetadata.ExtraEphys.Parameters.SampleRate;

    % xml making
    % Defaults... necessary at this point only for purposes of making an .xml l
%     SessionMetadata.ExtraEphys.Parameters.LfpSampleRate = LfpSampleRate;%from manual input
%     SessionMetadata.EphysDefaults.PointsPerWaveform = PointsPerWaveform;%from manual input
%     SessionMetadata.EphysDefaults.PeakPointinWaveform = PeakPointInWaveform;%from manual input
%     SessionMetadata.EphysDefaults.FeaturesPerWave = FeaturesPerWave;%from manual input
    %% Make XML for animal
    % aname = AnimalMetadata.AnimalName;
    % apath = AnimalMetadata.AnimalBasepath;
    pfiles = SessionMetadata.AnimalMetadata.Probes.ProbeLayoutFilenames;
    plugord = SessionMetadata.AnimalMetadata.Probes.PluggingOrder;
    bz_MakeXMLFromProbeMaps(basepath,basename,pfiles,plugord,SessionMetadata.ExtraEphys.Parameters);
end

%% Re-present data to the user to allow them to add new bad channels and bad shanks, after seeing the 


%% Save all metadata
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
