function bz_SessionNotesTemplate(basepath,AnimalMetadata)

if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

%%
if ~exist('AnimalMetadata','var')
    load(fullfile(basepath,[basename '.AnimalMetadata.mat']))
end
SessionNotes.basepath = basepath;
SessionNotes.basename = basename;

%% HUMAN INPUT BELOW

% Surgery and Animal metadata
if AnimalMetadata.Modules.AnimalAndSurgery
    SessionNotes.Animal.WeightGrams = [250];
end

% Extracell Ephys metadata
if AnimalMetadata.Modules.ExtracellEphys
    SessionNotes.ExtracellEphys.NumberOfTurnsSinceSurgery = [0 0];%vector, one entry per probe
    SessionNotes.ExtracellEphys.Probes.PluggingOrder     = [];%vector, one entry per probe. blank defaults to animal plugging order
    SessionNotes.ExtracellEphys.BadShanks = 'FromAnimalMetadata';%or [1 2 5] vector for this recording.  If value is "FromAnimalMetadata", field will be populated from AnimalMetadata.ExtracellEphys.CurrentBadChannels
         % These bad shanks will be used to populate bad channels
    SessionNotes.ExtracellEphys.BadChannels = 'FromAnimalMetadata';%or [1 2 5] vector for this recording.  If value is "FromAnimalMetadata", field will be populated from AnimalMetadata.ExtracellEphys.CurrentBadChannels
    SessionNotes.ExtracellEphys.ChannelNotes = {''};

    SessionNotes.ExtracellEphys.Parameters.LfpSampleRate = 1250;%assumed default
    SessionNotes.ExtracellEphys.Parameters.PointsPerWaveform = 32;%default
    SessionNotes.ExtracellEphys.Parameters.PeakPointInWaveform = 16;%default
    SessionNotes.ExtracellEphys.Parameters.FeaturesPerWave = 4;%default
end

%% Virus metadata
if AnimalMetadata.Modules.Virus
    %need to fill this in
end

%% Intracell ephys metadata
if AnimalMetadata.Modules.IntracellEphys
    %need to fill this in
end

%% Imaging metadata
if AnimalMetadata.Modules.Imaging
    %need to fill this in
end

%% Other for ?
SessionNotes.Other = {''};


%% Auto save
save(fullfile(basepath,[basename '.SessionNotes.mat']),'SessionNotes')

