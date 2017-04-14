function bz_SessionNotesTemplate(basepath,basename,AnimalMetadata)

if ~exist('AnimalMetadata','var')
    load(fullfile(basepath,[basename '_AnimalMetadata.mat']))
end

%%
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
    SessionNotes.ExtracellEphys.BadShanks = [];
    SessionNotes.ExtracellEphys.BadChannels = []; %0-indexed as in neuroscope
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
save(fullfile(basepath,[basename '_SessionNotes.mat']),'SessionNotes')

