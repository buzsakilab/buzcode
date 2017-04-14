function AnimalMetadata = bz_AnimalMetadataTextTemplate(basepath)
% This .m file is a generic template and will be copied to animal
% folders/session folders.  Edits here will change what is automatically
% copied to animal folders when bz_EditAnimalMetadta.m is called.
% 
% Developers should edit bz_AnimalMetadataTextTemplate.m to change how
% AnimalMetadata is created.
% 
% Brendon Watson 2017


%% Initial variable parsing
if ~exist('bedasepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

%% HUMAN INPUT SECTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name and path info
AnimalMetadata.AnimalName = 'Ket1';
AnimalMetadata.AnimalBasepath = '/balrog_zpool/Ket1';%this can be changed for each computer and can then act as a handle for all subsequent analyses

% Determine which modules to use for this animal... 1 for use, 0 for not use 
AnimalMetadata.Modules.AnimalAndSurgery = 1;
AnimalMetadata.Modules.ExtracellEphys = 1;
AnimalMetadata.Modules.Optogenetics = 0;
AnimalMetadata.Modules.Virus = 0;
AnimalMetadata.Modules.IntracellEphys = 0;
AnimalMetadata.Modules.Imaging = 0;

% Surgery and Animal metadata
if AnimalMetadata.Modules.AnimalAndSurgery
    AnimalMetadata.Animal.Species = 'Rat';
    AnimalMetadata.Animal.Strain = 'SpragueDawley';
    AnimalMetadata.Animal.GeneticLine = 'WildType';
    AnimalMetadata.Animal.Sex = 'Male';
    AnimalMetadata.Animal.DateOfBirth = '20161221';%YYYYMMDD format
    AnimalMetadata.Animal.WeightGramsAtSurgery = 405;%grams

    AnimalMetadata.Surgery.Date = '20170317';
    AnimalMetadata.Surgery.Anesthesia.Name = 'Isoflurane';
    AnimalMetadata.Surgery.Anesthesia.ConcentrationPercent = '1';
    AnimalMetadata.Surgery.Analgesic.Name = 'Buprenex';
    AnimalMetadata.Surgery.Analgesic.Milligrams = 0.06;%usually given at 0.15mg/ml
    AnimalMetadata.Surgery.Antibiotics.Topical = 'Neopredef';
    AnimalMetadata.Surgery.Antibiotics.Intraperitoneal = '';
    AnimalMetadata.Surgery.Complications = '';
    AnimalMetadata.Surgery.DamageSites = '';
    AnimalMetadata.Surgery.Notes = 'Good';
end

% Extracell Ephys metadata
if AnimalMetadata.Modules.ExtracellEphys
    %On probe subfields below: if multiple probes, put in one entry in each field
    %per probe, make sure they align with each other properly and all
    %subsequent assumptions will work.
    AnimalMetadata.ExtracellEphys.Probes.UmPerScrewTurn = [288 288];
    AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes = 2;
    AnimalMetadata.ExtracellEphys.Probes.TargetRegions = {'dCA1','mPFC'};
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Anteroposterior = [-3.5,2.7];%one for each probe
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Mediolateral = [2.5,0.3];
    AnimalMetadata.ExtracellEphys.Probes.ImplantAngle.Anteroposterior = [0,0];%degrees of top anterior as sitting behind animal
    AnimalMetadata.ExtracellEphys.Probes.ImplantAngle.Mediolateral = [0,10];%degrees clockwise as sitting behind animal
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.DepthFromSurface = [1.5,2];
    AnimalMetadata.ExtracellEphys.Probes.OrientationOfProbe.FirstGroupRelativeToLastGroupClockwiseDegreesAnteriorIsZero = [90,135];%assumes linear arrays
    AnimalMetadata.ExtracellEphys.Probes.OrientationOfProbe.GroupOffsetsFromCenter_ApMlDv = [];%for non-linear arrangements: group x 3 coordinates for change from center
    AnimalMetadata.ExtracellEphys.Probes.PluggingOrder = [1,2];% order will be represented in .xml, ie if intan splitter dicates
    AnimalMetadata.ExtracellEphys.Probes.SiteSizesInUmSq = [160];%In square microns
    AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames = {'NRX_Buzsaki64_5X12';'NRX_Buzsaki64_8X8'};%filenames in /buzcode/GeneralComputation/geometries
    AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilenames = {'Ket1_Impedances_5Shank.csv','Ket1_Impedances_8Shank.csv'};%Filenames in basepath folder, or leave as {} if none

%     AnimalMetadata.ExtracellEphys.CurrentBadChannels = [];
%     AnimalMetadata.ExtracellEphys.CurrentBadShanks = [];%not used now
    
    AnimalMetadata.ExtracellEphys.Parameters.SampleRate = 20000;%Lab default
    AnimalMetadata.ExtracellEphys.Parameters.Amplification = 1;%Intan digitized on chip, let's say 1
    AnimalMetadata.ExtracellEphys.Parameters.VoltsPerUnit = 0.0000002;%Intan default                
    AnimalMetadata.ExtracellEphys.Parameters.BitsPerSample = 16;%Intan default
    AnimalMetadata.ExtracellEphys.Parameters.VoltageRange = 10;%not used except to make xml
    AnimalMetadata.ExtracellEphys.Parameters.LfpSampleRate = 1250;%Usual desired default
    AnimalMetadata.ExtracellEphys.Parameters.PointsPerWaveform = 32;%default
    AnimalMetadata.ExtracellEphys.Parameters.PeakPointInWaveform = 16;%default
    AnimalMetadata.ExtracellEphys.Parameters.FeaturesPerWave = 4;%default

    % SessionMetadata.ExtracellEphys.Parameters.NumberOfChannels = 64;
    % This is actually set later in bz_SetSessionMetadata, based on ProbeLayoutFilenames
end

% Optogenetics metadata
if AnimalMetadata.Modules.Optogenetics
    AnimalMetadata.Optogenetics.Opsin = [];    
    FiberNum = 1;%copy more of these blocks - one per probe
    AnimalMetadata.Optogenetics.Fibers(ProbeNum).FiberRadiusMicrons = [];
    AnimalMetadata.Optogenetics.Fibers(ProbeNum).LightSourceWavelengthNm = [];
    AnimalMetadata.Optogenetics.Fibers(ProbeNum).LightSourceDevice = [];
    AnimalMetadata.Optogenetics.Fibers(FiberNum).AffiliatedProbeNumber = [];%Probe if proximal to a probe
    AnimalMetadata.Optogenetics.Fibers(FiberNum).AffiliatedShankSpikeGroupIndex = [];%If proximal to a shank 
    AnimalMetadata.Optogenetics.Fibers(FiberNum).DistanceProixmalToShankTipInUm = [];
    AnimalMetadata.Optogenetics.Fibers(FiberNum).APCoordinates = [];%if if not with a probe
    AnimalMetadata.Optogenetics.Fibers(FiberNum).MLCoordinates = [];
    AnimalMetadata.Optogenetics.Fibers(FiberNum).DVCoordinates = [];
%     AnimalMetadata.Optogenetics.FibersOnProbes(ProbeNum).NumFibersPerProbe = [];
%     AnimalMetadata.Optogenetics.FibersOnProbes(ProbeNum).FiberRadiusMicrons = [];
%     AnimalMetadata.Optogenetics.FibersOnProbes(ProbeNum).LightSourceWavelengthNm = [];
%     AnimalMetadata.Optogenetics.FibersOnProbes(ProbeNum).LightSourceDevice = [];
%     AnimalMetadata.Optogenetics.FibersOnProbes(ProbeNum).ShankNumbers = [];
%     AnimalMetadata.Optogenetics.FibersOnProbes(ProbeNum).DistanceProixmalToShankTip = [];
% 
%     AnimalMetadata.Optogenetics.FibersAwayFromProbes.FiberRadiusMicrons = [];
%     AnimalMetadata.Optogenetics.FibersAwayFromProbes.LightSourceWavelengthNm = [];
%     AnimalMetadata.Optogenetics.FibersAwayFromProbes.LightSourceDevice = [];
%     AnimalMetadata.Optogenetics.FibersAwayFromProbes.APCoordinates = [];
%     AnimalMetadata.Optogenetics.FibersAwayFromProbes.MLCoordinates = [];
%     AnimalMetadata.Optogenetics.FibersAwayFromProbes.DVCoordinates = [];
    
    %need to fill this in
end

% Virus metadata
if AnimalMetadata.Modules.Virus
    AnimalMetadata.Virus.Strain = '';
    AnimalMetadata.Virus.Coordinates.Anteroposterior = [];%one for each injection
    AnimalMetadata.Virus.InjectionDate = '';
end

% Intracell ephys metadata
if AnimalMetadata.Modules.IntracellEphys
    %need to fill this in
end

% Imaging metadata
if AnimalMetadata.Modules.Imaging
    %need to fill this in
end

% Other for ?
AnimalMetadata.Other = {''};

%% Automated after this point, depending on modules used
if AnimalMetadata.Modules.ExtracellEphys
    [PerGroupSuperficialToDeep,SpatialXY,NumChansPerProbe,GroupsPerChannel] = bz_ReadProbeMapFiles(AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames);
    AnimalMetadata.ExtracellEphys.Probes.NumGroupsPerProbe = sum(~cellfun(@isempty,PerGroupSuperficialToDeep),2);
    AnimalMetadata.ExtracellEphys.Probes.WithinProbeXYLocations = SpatialXY;
    AnimalMetadata.ExtracellEphys.Probes.NumChansPerProbe = NumChansPerProbe;%to do
    AnimalMetadata.ExtracellEphys.Probes.ProbeSpikeGroupLayoutSuperficialToDeep = PerGroupSuperficialToDeep;
    AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal = sum(NumChansPerProbe);
    %make lookup tables for probe number and anatomy for each channel
    po = AnimalMetadata.ExtracellEphys.Probes.PluggingOrder;
    lut = [];
    lut_ap = [];
    glut = [];
    glut_ap = [];
    pglut = [];
    pglut_p = [];
    for pidx = 1:AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes
        lut = cat(1,lut,pidx*ones(AnimalMetadata.ExtracellEphys.Probes.NumChansPerProbe(pidx),1));
        lut_ap = cat(1,lut_ap,po(pidx)*ones(AnimalMetadata.ExtracellEphys.Probes.NumChansPerProbe(po(pidx)),1));
        glut = cat(1,glut,GroupsPerChannel{pidx}+length(pglut));
        glut_ap = cat(1,glut_ap,GroupsPerChannel{po(pidx)});
        pglut = cat(1,pglut,[1:AnimalMetadata.ExtracellEphys.Probes.NumGroupsPerProbe(pidx)]'+length(pglut));
        pglut_p = cat(1,pglut_p,pidx*ones(AnimalMetadata.ExtracellEphys.Probes.NumGroupsPerProbe(pidx),1));
    end

    AnimalMetadata.ExtracellEphys.Probes.ProbeToGroupLookupTable.Labels = {'ProbeNumber';'GroupNumber'};
    AnimalMetadata.ExtracellEphys.Probes.ProbeToGroupLookupTable.Table = [pglut_p pglut];
    AnimalMetadata.ExtracellEphys.Channels.ChannelToProbeLookupTable.Labels = {'ChannelNumber';'ProbeNumber'};
    AnimalMetadata.ExtracellEphys.Channels.ChannelToProbeLookupTable.Table = [[1:AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal]' lut];
    AnimalMetadata.ExtracellEphys.Channels.ChannelToGroupLookupTable.Labels = {'ChannelNumber';'GroupNumber'};
    AnimalMetadata.ExtracellEphys.Channels.ChannelToGroupLookupTable.Table = [[1:AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal]' glut];
    AnimalMetadata.ExtracellEphys.Channels.ChannelToAnatomyLookupTable.Labels = {'AnatomyNameIndexedByChannelNumber'};
    AnimalMetadata.ExtracellEphys.Channels.ChannelToAnatomyLookupTable.Table = AnimalMetadata.ExtracellEphys.Probes.TargetRegions(lut)';
    AnimalMetadata.ExtracellEphys.Channels.ChannelToProbeLookupTable_AsPlugged.Labels = {'ChannelNumber';'ProbeNumber'};
    AnimalMetadata.ExtracellEphys.Channels.ChannelToProbeLookupTable_AsPlugged.Table = [[1:AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal]' lut_ap];
    AnimalMetadata.ExtracellEphys.Channels.ChannelToGroupLookupTable_AsPlugged.Labels = {'ChannelNumber';'GroupNumber'};
    AnimalMetadata.ExtracellEphys.Channels.ChannelToGroupLookupTable_AsPlugged.Table = [[1:AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal]' glut_ap];
    AnimalMetadata.ExtracellEphys.Channels.ChannelToAnatomyLookupTable_AsPlugged.Labels = {'AnatomyNameIndexedByChannelNumber'};
    AnimalMetadata.ExtracellEphys.Channels.ChannelToAnatomyLookupTable_AsPlugged.Table = AnimalMetadata.ExtracellEphys.Probes.TargetRegions(lut_ap)';

    %Get impedances per channel based on intan impedance files
    if ~isempty(AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilenames)
        AnimalMetadata.ExtracellEphys.Channels.ImpedanceByChannel = cell(AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes);
        for pidx = 1:length(AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilenames)
            tf = fullfile(basepath,AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilenames{pidx});
            txt = read_mixed_csv(tf);
            numchans = size(txt,1)-1;
        %     channums = 1:numchans;
            for cidx = 1:numchans
                eloc = strfind(txt{cidx+1,5},'e');
                n1 = str2num(txt{cidx+1,5}(1:eloc-1));
                n2 = str2num(txt{cidx+1,5}(eloc+1:end));
                imp(cidx) = n1*10^n2;
            end
        %     AnimalMetadata.ExtracellEphys.Probes.ImpedanceByChannel{pidx} = cat(2,channums',imp');
            AnimalMetadata.ExtracellEphys.Channels.ImpedanceByChannel{pidx} = imp';
        end
    end

    % someone should combine XY with probe orientation angle and implant
   %coordinates to map each site at day of implant

    % fix defaults
    AnimalMetadata.EphysDefaults.NumberOfChannels = AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal;

end

%% Save
save(fullfile(basepath,[basename '.AnimalMetadata.mat']),'AnimalMetadata')


function lineArray = read_mixed_csv(fileName,delimiter)
% copied from http://stackoverflow.com/questions/4747834/import-csv-file-with-mixed-data-types
% Brendon Watson 2014

if ~exist('delimiter','var')
    delimiter = ',';
end

  fid = fopen(fileName,'r');   %# Open the file
  lineArray = cell(100,1);     %# Preallocate a cell array (ideally slightly
                               %#   larger than is needed)
  lineIndex = 1;               %# Index of cell to place the next line in
  nextLine = fgetl(fid);       %# Read the first line from the file
  while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
    lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
    lineIndex = lineIndex+1;          %# Increment the line index
    nextLine = fgetl(fid);            %# Read the next line from the file
  end
  fclose(fid);                 %# Close the file
  lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
  for iLine = 1:lineIndex-1              %# Loop over lines
    lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                        'Delimiter',delimiter);
    lineData = lineData{1};              %# Remove cell encapsulation
    if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
      lineData{end+1} = '';                     %#   ends with a delimiter
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
  end