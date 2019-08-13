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
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

%% HUMAN INPUT SECTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name and path info
AnimalMetadata.AnimalName = 'Jasper';
AnimalMetadata.AnimalBasepath = 'F:\Jasper';%this can be changed for each computer and can then act as a handle for all subsequent analyses

% Determine which modules to use for this animal... 1 for use, 0 for not use 
AnimalMetadata.Modules.AnimalAndSurgery = true;
AnimalMetadata.Modules.ExtracellEphys = true;
AnimalMetadata.Modules.Optogenetics = false;
AnimalMetadata.Modules.Virus = false;
AnimalMetadata.Modules.IntracellEphys = false;
AnimalMetadata.Modules.Imaging = false;

% Surgery and Animal metadata
if AnimalMetadata.Modules.AnimalAndSurgery
    AnimalMetadata.Animal.Species = 'Rat';
    AnimalMetadata.Animal.Strain = 'SpragueDawley';
    AnimalMetadata.Animal.GeneticLine = 'WildType';
    AnimalMetadata.Animal.Sex = 'Male';
    AnimalMetadata.Animal.DateOfBirth = '20181112';%YYYYMMDD format
    AnimalMetadata.Animal.WeightGramsAtSurgery = 350;%grams

    AnimalMetadata.Surgery.Date = '20190221';%YYYYMMDD format
    AnimalMetadata.Surgery.Anesthesia.Name = 'Isoflurane';
    AnimalMetadata.Surgery.Anesthesia.ConcentrationPercent = '';
    
    AnimalMetadata.Surgery.Analgesic.Name = 'Carprofen';
    AnimalMetadata.Surgery.Analgesic.DoseMgPerKg = [5];
    AnimalMetadata.Surgery.Analgesic.Milligrams = AnimalMetadata.Surgery.Analgesic.DoseMgPerKg/1000*AnimalMetadata.Animal.WeightGramsAtSurgery;%usually given at 0.15mg/ml
    AnimalMetadata.Surgery.Steroids.Name = 'Methylprednisolone';
    AnimalMetadata.Surgery.Steroids.DoseMgPerKg = [30];
    AnimalMetadata.Surgery.Steroids.Name = AnimalMetadata.Surgery.Steroids.DoseMgPerKg /1000 * AnimalMetadata.Animal.WeightGramsAtSurgery;
    AnimalMetadata.Surgery.Antibiotics.Intraperitoneal.Name = '';
    AnimalMetadata.Surgery.Antibiotics.Intraperitoneal.DoseMgPerKg = [];
    AnimalMetadata.Surgery.Antibiotics.Intraperitoneal.Milligrams = AnimalMetadata.Surgery.Antibiotics.Intraperitoneal.DoseMgPerKg /1000 * AnimalMetadata.Animal.WeightGramsAtSurgery;;
    AnimalMetadata.Surgery.Antibiotics.Subcutaneous.Name = 'Cefazolin';
    AnimalMetadata.Surgery.Antibiotics.Subcutaneous.DoseMgPerKg = [60];
    AnimalMetadata.Surgery.Antibiotics.Subcutaneous.Milligrams = AnimalMetadata.Surgery.Antibiotics.Subcutaneous.DoseMgPerKg /1000 * AnimalMetadata.Animal.WeightGramsAtSurgery;;
    AnimalMetadata.Surgery.Antibiotics.SkinTopical.Name = '';
    
    AnimalMetadata.Surgery.DamageSites = '';
    AnimalMetadata.Surgery.Complications = '';
    AnimalMetadata.Surgery.Notes = '';
end

% Extracell Ephys metadata
if AnimalMetadata.Modules.ExtracellEphys
    %On probe subfields below: if multiple probes, put in one entry in each field
    %per probe, make sure they align with each other properly and all
    %subsequent assumptions will work.
    AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes = 2;
    AnimalMetadata.ExtracellEphys.Probes.TargetRegions = {'mPFC','dCA1'};
    AnimalMetadata.ExtracellEphys.Probes.TargetHemisphere = {'right','right'};
    AnimalMetadata.ExtracellEphys.Probes.UmPerScrewTurn = [205 205];
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Anteroposterior = [3.5,-3.8];%one for each probe
    AnimalMetadata.ExtracellEphys.Probes.ImplantAngle.Anteroposterior = [0,0];%degrees of top anterior as sitting behind animal
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Mediolateral = [-0.6,-3,5];
    AnimalMetadata.ExtracellEphys.Probes.ImplantAngle.Mediolateral = [0,0];%degrees clockwise as sitting behind animal
    AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.DepthFromSurface = [-2.5,-1.8];
    AnimalMetadata.ExtracellEphys.Probes.OrientationOfProbe.FirstGroupRelativeToLastGroupClockwiseDegreesAnteriorIsZero = [270,270];%assumes linear arrays
    AnimalMetadata.ExtracellEphys.Probes.OrientationOfProbe.GroupOffsetsFromCenter_ApMlDv = [];%for non-linear arrangements: group x 3 coordinates for change from center
    AnimalMetadata.ExtracellEphys.Probes.PluggingOrder = [1,2];% order will be represented in .xml, ie if intan splitter dicates
    AnimalMetadata.ExtracellEphys.Probes.SiteSizesInUmSq = [177,177];%In square microns
    AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames = {'NRX_Buzsaki64_8X8';'NRX_Buzsaki64_5X12'};%filenames in /buzcode/GeneralComputation/geometries
    %ONLY ENTER THIS MANUALLY IF YOU WOULD LIKE TO USE THE EXISTING .xml 
    %SPIKE GROUPS INSTEAD OF CALCULATING FROM A PROBE GEOMETRY FILE
    %Enter a list of spike group numbers (from neuroscope) for each probe
    AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups = {[],[]};

%     AnimalMetadata.ExtracellEphys.Probes.BadShanks = [];%not used now
    AnimalMetadata.ExtracellEphys.Channels.BadChannels = [];

    AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilepaths = {};%Filenames in basepath folder, or leave as {} if none
    
    %For setting up to record extra channels for opto, frames etc.  
    % NumExtraChansPerExtraGroup field will be a vector of number of channels per
    % extra group.   The number of groups is set by the number of entries... 
    % the number within each group is signified by the number in each
    % entry.  (ie [3 5] signifies 2 post-probe groups, first with 3 channels, the 
    % second with 5 channels.
    AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup = [0];
    AnimalMetadata.ExtracellEphys.ExtraChannels.GroupNames = {};%one per group
    AnimalMetadata.ExtracellEphys.ExtraChannels.ChannelNames = {};%one per channel

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
    AnimalMetadata.ExtracellEphys.Parameters.FeaturesPerWave = 3;%for PCA - default

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
    GroupsFromXml = false;
    GroupsFromXml = ~isempty(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups) && ...
        ~isempty(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups{1}) && ...
        isnumeric(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikegroups{1});
    if GroupsFromXml
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

    else     %If you want to build the spikegroups from the geometry files
        display('Reading probe layout info from the probe map file...')
        [PerGroupSuperficialToDeep,SpatialXY,NumChansPerProbe,GroupsPerChannel] = bz_ReadProbeGeometryFiles...
            (AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames,AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup);
        AnimalMetadata.ExtracellEphys.Probes.WithinProbeXYLocations = SpatialXY;
    end

    %handling group names - for case of extra channels - esp if not named
    %by user
    for gnidx = 1:length(AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup)
        addsubstitutename = false;
        if length(AnimalMetadata.ExtracellEphys.ExtraChannels.GroupNames) >= gnidx
            if isempty(AnimalMetadata.ExtracellEphys.ExtraChannels.GroupNames)
                addsubstitutename = true;
            else
            end
        else 
            addsubstitutename = true;
        end
        if addsubstitutename
            AnimalMetadata.ExtracellEphys.ExtraChannels.GroupNames{gnidx} = ['ExtraChannelGroup' num2str(gnidx)];
        end
    end
    AnimalMetadata.ExtracellEphys.Probes.TargetRegions = cat(2,...
        AnimalMetadata.ExtracellEphys.Probes.TargetRegions,...
        AnimalMetadata.ExtracellEphys.ExtraChannels.GroupNames);


    
    AnimalMetadata.ExtracellEphys.Probes.NumGroupsPerProbe = sum(~cellfun(@isempty,PerGroupSuperficialToDeep),2);
    AnimalMetadata.ExtracellEphys.Probes.NumChansPerProbe = NumChansPerProbe;%to do
    AnimalMetadata.ExtracellEphys.Probes.ProbeSpikeGroupLayoutSuperficialToDeep = PerGroupSuperficialToDeep;
    AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal = sum(NumChansPerProbe);
    
    % sessionInfo Style Groups
    SpikeGroups.nGroups = sum(AnimalMetadata.ExtracellEphys.Probes.NumGroupsPerProbe);
    SpikeGroups.nSamples = AnimalMetadata.ExtracellEphys.Parameters.PointsPerWaveform * ones(1,SpikeGroups.nGroups);
    SpikeGroups.groups = {};
    for pidx = 1:size(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikeGroupLayoutSuperficialToDeep,1)       
        for gidx = 1:size(AnimalMetadata.ExtracellEphys.Probes.ProbeSpikeGroupLayoutSuperficialToDeep,2)
            tgroup = AnimalMetadata.ExtracellEphys.Probes.ProbeSpikeGroupLayoutSuperficialToDeep{pidx,gidx};
            if ~isempty(tgroup)
                if prod(size(tgroup))>0
                    SpikeGroups.groups = cat(2,SpikeGroups.groups,tgroup);
                    SpikeGroups.groups{end} = SpikeGroups.groups{end}';
                end
            end
        end
    end
    AnimalMetadata.ExtracellEphys.SpikeGroups = SpikeGroups;
    
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
        glut = cat(1,glut,GroupsPerChannel{pidx});
        glut_ap = cat(1,glut_ap,GroupsPerChannel{po(pidx)});
        pglut = cat(1,pglut,[1:AnimalMetadata.ExtracellEphys.Probes.NumGroupsPerProbe(pidx)]'+length(pglut));
        pglut_p = cat(1,pglut_p,pidx*ones(AnimalMetadata.ExtracellEphys.Probes.NumGroupsPerProbe(pidx),1));
    end
    numrealprobes = AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes;
    numspikegroups = AnimalMetadata.ExtracellEphys.SpikeGroups.nGroups;
    if ~isempty(AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup);
       for eidx = 1:length(AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes);
          lut = cat(1,lut,(eidx+numrealprobes)*ones(AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup(eidx),1));
          lut_ap = cat(1,lut_ap,(eidx+numrealprobes)*ones(AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup(eidx),1));
          glut = cat(1,glut,(eidx+max(glut))*ones(AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup(eidx),1));
          glut_ap = cat(1,glut_ap,(eidx+max(glut_ap))*ones(AnimalMetadata.ExtracellEphys.ExtraChannels.NumExtraChansPerExtraGroup(eidx),1));
          
          pglut = cat(1,pglut,1+length(pglut));%one "probe" for each set of extra channels
          pglut_p = cat(1,pglut_p,(eidx+numrealprobes));
       end
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
    if ~isempty(AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilepaths)
        AnimalMetadata.ExtracellEphys.Channels.ImpedanceByChannel = cell(AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes);
        for pidx = 1:length(AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilepaths)
            tf = fullfile(basepath,AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilepaths{pidx});
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

   % someone should combine XY with probe orientation angle and implant...
   % to get full stereotaxic coordinates for each site... wouldn't that be
   % nice
   
    AnimalMetadata.ExtracellEphys.NumberOfChannels = AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal;
end


%% Make session info for ease (3/2019)
%should change this to use bz_editSessionInfo and use these as inputs to
%that
%maybe even abstract wrapper so field name and field value go into
%bz_editSession info as name-value pairs?
sessionInfo.session.name = basename;
sessionInfo.session.path = basepath;
sessionInfo.spikeGroups = AnimalMetadata.ExtracellEphys.SpikeGroups;
sessionInfo.nChannels = AnimalMetadata.ExtracellEphys.NumberOfChannels;
sessionInfo.channels = [0:AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal-1];
sessionInfo.nBits = AnimalMetadata.ExtracellEphys.Parameters.BitsPerSample;
sessionInfo.rates.lfp = AnimalMetadata.ExtracellEphys.Parameters.LfpSampleRate;
sessionInfo.rates.wideband = AnimalMetadata.ExtracellEphys.Parameters.SampleRate;
sessionInfo.rates.video = 0;%default
sessionInfo.FileName = basename;
sessionInfo.SampleTime = 1/AnimalMetadata.ExtracellEphys.Parameters.SampleRate* 1e+6;%us apparently
sessionInfo.nElecGps = [];
for gidx = 1:(sessionInfo.spikeGroups.nGroups)%making elecgp... this is dumb
    tchans = sessionInfo.spikeGroups.groups{gidx};
    for cidx = 1:length(tchans)
        sessionInfo.ElecGp{gidx}.channel{cidx} = num2str(tchans(cidx));
    end
end
sessionInfo.HiPassFreq = [];
sessionInfo.lfpSampleRate = AnimalMetadata.ExtracellEphys.Parameters.LfpSampleRate;
sessionInfo.VoltageRange = AnimalMetadata.ExtracellEphys.Parameters.VoltageRange;
sessionInfo.Amplification = AnimalMetadata.ExtracellEphys.Parameters.Amplification;
sessionInfo.Offset = nan;
d = AnimalMetadata.Surgery.Date;
sessionInfo.Date = strcat(d(1:4),'-',d(5:6),'-',d(7:8));
for gidx = 1:(sessionInfo.spikeGroups.nGroups)%making AnatGrps
    sessionInfo.AnatGrps(gidx).Channels = sessionInfo.spikeGroups.groups{gidx};
end
for gidx = 1:(sessionInfo.spikeGroups.nGroups)%making SpkGrps
    sessionInfo.SpkGrps(gidx).Channels = sessionInfo.spikeGroups.groups{gidx};
    sessionInfo.SpkGrps(gidx).nSamples = AnimalMetadata.ExtracellEphys.Parameters.PointsPerWaveform;
    sessionInfo.SpkGrps(gidx).peakSample = AnimalMetadata.ExtracellEphys.Parameters.PeakPointInWaveform;
    sessionInfo.SpkGrps(gidx).nFeatures = AnimalMetadata.ExtracellEphys.Parameters.FeaturesPerWave;
end
%sessionInfo.Units = 
sessionInfo.region = AnimalMetadata.ExtracellEphys.Channels.ChannelToAnatomyLookupTable.Table';
sessionInfo.badchannels = AnimalMetadata.ExtracellEphys.Channels.BadChannels;
sessionInfo.badshanks = AnimalMetadata.ExtracellEphys.CurrentBadShanks;%not used now

AnimalMetadata.sessionInfo = sessionInfo;

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