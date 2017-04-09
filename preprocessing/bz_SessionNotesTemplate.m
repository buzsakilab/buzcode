function bz_SessionNotesTemplate(basepath,basename,AnimalMetadata)

SessionNotes.basepath = basepath;
SessionNotes.basename = basename;

%% HUMAN INPUT BELOW

% Surgery and Animal metadata
if AnimalMetadata.Modules.AnimalAndSurgery
    SessionNotes.Animal.WeightGrams = [250];
end

% Extracell Ephys metadata
if AnimalMetadata.Modules.ExtracellEphys
    SessionNotes.ExtraEphys.NumberOfTurnsSinceSurgery = [0 0];%vector, one entry per probe
    SessionNotes.ExtraEphys.BadShanks = [];
    SessionNotes.ExtraEphys.BadChannels = [];
    SessionNotes.ExtraEphys.ChannelNotes = {''};
    
    SessionNotes.ExtraEphys.Parameters.LfpSampleRate = 1250;%assumed default
    SessionNotes.ExtraEphys.Parameters.PointsPerWaveform = 32;%default
    SessionNotes.ExtraEphys.Parameters.PeakPointInWaveform = 16;%default
    SessionNotes.ExtraEphys.Parameters.FeaturesPerWave = 4;%default
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

