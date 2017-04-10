function bz_AnimalNotesTemplate(basepath,basename)

%% Variable parsing
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

%% HUMAN INPUT BELOW
% Name and path info
AnimalNotes.AnimalName = 'Ket1';
AnimalNotes.AnimalBasepath = '/balrog_zpool/Ket1';%this can be changed for each computer and can then act as a handle for all subsequent analyses

% Determine which modules to use for this animal... 1 for use, 0 for not use 
AnimalNotes.Modules.AnimalAndSurgery = 1;
AnimalNotes.Modules.ExtracellEphys = 1;
AnimalNotes.Modules.Virus = 0;
AnimalNotes.Modules.IntracellEphys = 0;
AnimalNotes.Modules.Imaging = 0;

% Surgery and Animal metadata
if AnimalNotes.Modules.AnimalAndSurgery

    AnimalNotes.Animal.Species = 'Rat';
    AnimalNotes.Animal.Strain = 'SpragueDawley';
    AnimalNotes.Animal.GeneticLine = 'WildType';
    AnimalNotes.Animal.Sex = 'Male';
    AnimalNotes.Animal.DateOfBirth = '20161221';%YYYYMMDD format
    AnimalNotes.Animal.WeightGramsAtSurgery = 405;%grams

    AnimalNotes.Surgery.Date = '20170317';
    AnimalNotes.Surgery.Anesthesia.Name = 'Isoflurane';
    AnimalNotes.Surgery.Anesthesia.ConcentrationPercent = '1';
    AnimalNotes.Surgery.Analgesic.Name = 'Buprenex';
    AnimalNotes.Surgery.Analgesic.Milligrams = 0.06;%usually given at 0.15mg/ml
    AnimalNotes.Surgery.Antibiotics.Topical = 'Neopredef';
    AnimalNotes.Surgery.Antibiotics.Intraperitoneal = '';
    AnimalNotes.Surgery.Complications = '';
    AnimalNotes.Surgery.DamageSites = '';
    AnimalNotes.Surgery.Notes = 'Good';
end

% Extracell Ephys metadata
if AnimalNotes.Modules.ExtracellEphys
    %On probe subfields below: if multiple probes, put in one entry in each field
    %per probe, make sure they align with each other properly and all
    %subsequent assumptions will work.
    AnimalNotes.Probes.UmPerScrewTurn = [288 288];
    AnimalNotes.Probes.NumberOfProbes = 2;
    AnimalNotes.Probes.TargetRegions = {'dCA1','mPFC'};
    AnimalNotes.Probes.ImplantCoordinates.Anteroposterior = [-3.5,2.7];%one for each probe
    AnimalNotes.Probes.ImplantCoordinates.Mediolateral = [2.5,0.3];
    AnimalNotes.Probes.ImplantAngle.Anteroposterior = [0,0];%degrees of top anterior as sitting behind animal
    AnimalNotes.Probes.ImplantAngle.Mediolateral = [0,10];%degrees clockwise as sitting behind animal
    AnimalNotes.Probes.ImplantCoordinates.DepthFromSurface = [1.5,2];
    AnimalNotes.Probes.OrientationOfProbe.FirstGroupRelativeToLastGroupClockwiseDegreesAnteriorIsZero = [90,135];%assumes linear arrays
    AnimalNotes.Probes.OrientationOfProbe.GroupOffsetsFromCenter_ApMlDv = [];%for non-linear arrangements: group x 3 coordinates for change from center
    AnimalNotes.Probes.PluggingOrder = [2,1];% order will be represented in .xml, ie if intan splitter dicates
    AnimalNotes.Probes.SiteSizesInUmSq = [160];%In square microns
    AnimalNotes.Probes.ProbeLayoutFilenames = {'NRX_Buzsaki64_5X12';'NRX_Buzsaki64_8X8'};%filenames in /buzcode/GeneralComputation/geometries
    AnimalNotes.Channels.ImpedanceFilenames = {'Ket1_Impedances_5Shank.csv','Ket1_Impedances_8Shank.csv'};%Filenames in basepath folder, or leave as {} if none
end

%% Virus metadata
if AnimalNotes.Modules.Virus
    AnimalNotes.Virus.Strain = '';
    AnimalNotes.Virus.Coordinates.Anteroposterior = [];%one for each injection
    AnimalNotes.Virus.InjectionDate = '';
end

%% Intracell ephys metadata
if AnimalNotes.Modules.IntracellEphys
    %need to fill this in
end

%% Imaging metadata
if AnimalNotes.Modules.Imaging
    %need to fill this in
end

%% Other for ?
AnimalNotes.Other = {''};


%% Auto save
save(fullfile(basepath,[basename '_AnimalNotes.mat']),'AnimalNotes')

