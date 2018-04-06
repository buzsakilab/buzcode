function AnimalMetadata = bz_AnimalMetadataTextTemplate(basepath)
% This .m file is a generic template and will be copied to animal
% folders/session folders and will be called basename_AnimalMetadataText.m.
% It is represents code that will be used to
% create a basename.AnimalMetadata.mat file.  Edits here will be reflected
% in the .AnimalMetadata.mat created after this file is run.  
% 
% Developers interested in changing how AnimalMetadata is made generally 
% should edit bz_AnimalMetadataTextTemplate.m 
%
% INPUTS
%   basepath - computer path to folder where animal metadata is to go.
%   Often this is a path to an folder full of sessions all from one animal.
%
% OUTPUTS
%    AnimalMetadata and saved basename.AnimalMetadata.mat
%       - Struct array containing fields including surgical info, animal
%       name, viruses injected, probes used etc.
%
% See also: bz_EditAnimalMetadata, bz_RunAnimalMetadata, bz_SessionMetadataTextTemplate
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
AnimalMetadata.AnimalName = 'Mouse5';
AnimalMetadata.AnimalBasepath = 'C:\IntanRecording\MD\Mouse5\Mouse5_180130_150939';%this can be changed for each computer and can then act as a handle for all subsequent analyses

% Determine which modules to use for this animal... 1 for use, 0 for not use 
AnimalMetadata.Modules.AnimalAndSurgery = 1;
AnimalMetadata.Modules.ExtracellEphys = 1;
AnimalMetadata.Modules.Optogenetics = 0;
AnimalMetadata.Modules.Virus = 0;
AnimalMetadata.Modules.IntracellEphys = 0;
AnimalMetadata.Modules.Imaging = 0;

% Surgery and Animal metadata
if AnimalMetadata.Modules.AnimalAndSurgery
    AnimalMetadata.Animal.Species = 'Mouse';
    AnimalMetadata.Animal.Strain = 'C57Bl6';
    AnimalMetadata.Animal.GeneticLine = 'WildType';
    AnimalMetadata.Animal.Sex = 'Male';
    AnimalMetadata.Animal.DateOfBirth = '20171118';%YYYYMMDD format
    AnimalMetadata.Animal.WeightGramsAtSurgery = nan;%grams

    AnimalMetadata.Surgery.Date = '20180130';
    AnimalMetadata.Surgery.Anesthesia.Name = 'Isoflurane';
    AnimalMetadata.Surgery.Anesthesia.ConcentrationPercent = '';
    AnimalMetadata.Surgery.Analgesic.Name = '';
    AnimalMetadata.Surgery.Analgesic.Milligrams = nan;%usually given at 0.15mg/ml
    AnimalMetadata.Surgery.Antibiotics.Topical = '';
    AnimalMetadata.Surgery.Antibiotics.Intraperitoneal = '';
    AnimalMetadata.Surgery.Complications = '';
    AnimalMetadata.Surgery.DamageSites = '';
    AnimalMetadata.Surgery.Notes = '';
end

% Extracell Ephys metadata
if AnimalMetadata.Modules.ExtracellEphys
    %On probe subfields below: if multiple probes, put in one entry in each field
    %per probe, make sure they align with each other properly and all
    %subsequent assumptions will work.
    AnimalMetadata.ExtracellEphys.Probes.UmPerScrewTurn = [nan];
    AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes = 1;
    AnimalMetadata.ExtracellEphys.Probes.TargetRegions = {'S1'};
    AnimalMetadata.ExtracellEphys.Probes.TargetHemisphere = {'left'};
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Anteroposterior = [-1.5];%one for each probe
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Mediolateral = [3];
    AnimalMetadata.ExtracellEphys.Probes.ImplantAngle.Anteroposterior = [0];%degrees of top anterior as sitting behind animal
    AnimalMetadata.ExtracellEphys.Probes.ImplantAngle.Mediolateral = [-25];%degrees clockwise as sitting behind animal
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.DepthFromSurface = [1.1];
    AnimalMetadata.ExtracellEphys.Probes.OrientationOfProbe.FirstGroupRelativeToLastGroupClockwiseDegreesAnteriorIsZero = [nan];%assumes linear arrays
    AnimalMetadata.ExtracellEphys.Probes.OrientationOfProbe.GroupOffsetsFromCenter_ApMlDv = [];%for non-linear arrangements: group x 3 coordinates for change from center
    AnimalMetadata.ExtracellEphys.Probes.PluggingOrder = [1,2];% order will be represented in .xml, ie if intan splitter dicates
    AnimalMetadata.ExtracellEphys.Probes.SiteSizesInUmSq = [nan];%In square microns
    AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames = {'CED_H3_1X64'};%filenames in /buzcode/GeneralComputation/geometries
    AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilenames = {};%Filenames in basepath folder, or leave as {} if none

    %ONLY ENTER THIS MANUALLY IF YOU WOULD LIKE TO USE THE EXISTING .xml 
    %SPIKE GROUPS INSTEAD OF CALCULATING FROM A PROBE GEOMETRY FILE
    %Enter a list of spike group numbers (from neuroscope) for each probe
    AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups = {['spkgrps on probe1']};
    
    
%     AnimalMetadata.ExtracellEphys.CurrentBadChannels = [];
%     AnimalMetadata.ExtracellEphys.CurrentBadShanks = [];%not used now
    
    % Default ephys parameters for this anaimal.  These will be checked
    % against what was used on a per-session basis and warnings will be
    % issued if differences found.
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
    AnimalMetadata.Optogenetics.Opsin = {''};%one character string per probe   
    AnimalMetadata.Optogenetics.Promoter = {''};%one character string per probe
%copy more of these blocks - one per fiber
    FiberNum = 1;%copy more of these blocks - one per fiber
    AnimalMetadata.Optogenetics.Fibers(ProbeNum).FiberRadiusMicrons = [];
    AnimalMetadata.Optogenetics.Fibers(ProbeNum).LightSourceWavelengthNm = [];
    AnimalMetadata.Optogenetics.Fibers(ProbeNum).LightSourceDevice = [];
    AnimalMetadata.Optogenetics.Fibers(FiberNum).AffiliatedProbeNumber = [];%Probe if proximal to a probe. Probe number based on ordering above
    AnimalMetadata.Optogenetics.Fibers(FiberNum).AffiliatedShankSpikeGroupIndex = [];%If proximal to a shank 
    AnimalMetadata.Optogenetics.Fibers(FiberNum).DistanceProximalToShankTipInUm = [];
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
AnimalMetadata.ExtracellEphys.Probes.WithinProbeXYLocations = [];%declaring default since may not be declared below
if AnimalMetadata.Modules.ExtracellEphys
    %If you want to use .xml for spikegroups
    if isnumeric(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups{1})
        display('Getting spikegroups from the .xml...')
        %Load the XML for spike groups
        xmlname = fullfile(basepath,[basename,'.xml']);
        if ~exist(xmlname,'file')
            [FileName,PathName] = uigetfile('.xml','Find the .xml to load SpikeGroups from',basepath);
            xmlname = fullfile(PathName,FileName);
        end
        xmlparms = LoadParameters(xmlname);
        %Put Channels from each group into the PerGroupSuperficialToDeep matrix
        for pp = 1:length(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups)
            for gg = 1:length(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups{pp})
                thisgroup = AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups{pp}(gg);
                PerGroupSuperficialToDeep{pp,gg} = xmlparms.SpkGrps(thisgroup).Channels;
                GroupsPerChannel{pp}(PerGroupSuperficialToDeep{pp,gg}+1,1) = gg;
                GroupsPerChannel{pp}(GroupsPerChannel{pp}==0)=[];
            end
            NumChansPerProbe(pp) = length(cat(2,PerGroupSuperficialToDeep{pp,:}));
        end

    %If you want to build the spikegroups from the geometry
    elseif strcmp(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups{1},'spkgrps on probe1')
        display('Reading probe layout info from the probe map file...')
        [PerGroupSuperficialToDeep,SpatialXY,NumChansPerProbe,GroupsPerChannel] = bz_ReadProbeMapFiles(AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames);
        AnimalMetadata.ExtracellEphys.Probes.WithinProbeXYLocations = SpatialXY;
    else display('Spikegroup error')
    end
    
    AnimalMetadata.ExtracellEphys.Probes.NumGroupsPerProbe = sum(~cellfun(@isempty,PerGroupSuperficialToDeep),2);
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