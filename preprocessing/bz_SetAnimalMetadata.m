function AnimalMetaData = bz_SetAnimalMetadata(basepath,basename)
%MakeXMLFromProbeMaps - Generate a .xml file to accompany .dat files for a
%recording in the neuroscope/klusters/ndmanager system.  Uses a library of
%probe map layouts.
%
%  USAGE
%
%    MakeXMLFromProbeMaps(basepath,basename,ProbeFileName1,ProbeFileName2...)
%       Example:
%    MakeXMLFromProbeMaps(cd,'','NRX_Buzsaki64_8X8','NRX_Buzsaki64_6X10');
%
%    Writes a standardized .xml file based on a user-selection of probe
%    maps and in a sequence specified by the user (ie 64site probe first
%    then 32site probe second).  Probe maps can be found at:
%    /buzcode/tree/master/generalComputation/geometries
%
%  INPUT
%
%    basepath       Path to directory to which to write output xml file and
%                   where to potentially find .rhd file. Default is path to
%                   the current directory.
%
%    basename       Shared name for this file and all others for this
%                   recording.  Default will the name of the current
%                   folder.
%
%    ProbeFileNames Names of .xlsx files specifying probe geometries, must 
%                   be on the path, ie from
%                   buzcode/GeneralComputation/Geometries).
%                   Must have the following format features:
%                       - A column with row 1 having text "BY VERTICAL POSITION/SHANK (IE FOR DISPLAY)".
%                       In this row must be cells specifying SHANK 1, 
%                       SHANK 2, etc.  
%                       - Two rows right of that a column with channel
%                       numbers listed from superficial to deep in the
%                       shank group specified two to the left.
%
%  OUTPUT
%
%    (.xml file written to disk at basepath)
%
%  SEE
%
%    See also 
% Copyright (C) 2017 Brendon Watson
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

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

%% Value setting - humans must do this
AnimalMetaData.AnimalName = 'Ket1';
AnimalMetaData.AnimalBaseFolderPath = '/Users/brendon/Dropbox/Data/Ket1';%this can be changed for each computer and can then act as a handle for all subsequent analyses

AnimalMetaData.Species = 'Rat';
AnimalMetaData.Strain = 'SpragueDawley';
AnimalMetaData.GeneticLine = 'WildType';
AnimalMetaData.Sex = 'male';
AnimalMetaData.DateOfBirth = '20161221';%YYYYMMDD format
AnimalMetaData.WeightGramsAtSurgery = 405;%grams

AnimalMetaData.Surgery.Date = '20140317';
AnimalMetaData.Surgery.Anesthesia.Name = 'Isoflurane';
AnimalMetaData.Surgery.Anesthesia.ConcentrationPercent = '1';
AnimalMetaData.Surgery.Analgesic.Name = 'Buprenex';
AnimalMetaData.Surgery.Analgesic.Milligrams = 0.06;%usually given at 0.15mg/ml
AnimalMetaData.Surgery.Antibiotics.Topical = 'Neopredef';
AnimalMetaData.Surgery.Antibiotics.Intraperitoneal = '';
AnimalMetaData.Surgery.Complications = '';
AnimalMetaData.Surgery.DamageSites = '';
AnimalMetaData.Surgery.Notes = 'Good';

AnimalMetaData.Virus.Strain = '';
AnimalMetaData.Virus.Coordinates.Anteroposterior = [];%one for each injection
AnimalMetaData.Virus.InjectionDate = '';

AnimalMetaData.Probes.NumberOfProbes = 2;
AnimalMetaData.Probes.TargetAnatomy = {'dCA1','mPFC'};
AnimalMetaData.Probes.ImplantCoordinates.Anteroposterior = [-3.5,2.7];%one for each probe
AnimalMetaData.Probes.ImplantCoordinates.Mediolateral = [2.5,0.3];
AnimalMetaData.Probes.ImplantAngle.Anteroposterior = [0,0];%degrees of top anterior as sitting behind animal
AnimalMetaData.Probes.ImplantAngle.Mediolateral = [0,10];%degrees clockwise as sitting behind animal
AnimalMetaData.Probes.ImplantCoordinates.DepthFromSurface = [1.5,2];
AnimalMetaData.Probes.OrientationOfProbe.FirstGroupRelativeToLastGroupClockwiseDegrees = [90,135];%assumes linear arrays
AnimalMetaData.Probes.OrientationOfProbe.GroupOffsetsFromCenter_ApMlDv = [];%for non-linear arrangements: group x 3 coordinates for change from center
AnimalMetaData.Probes.PluggingOrder = [2,1];% order will be represented in .xml, ie if intan splitter dicates
AnimalMetaData.Probes.SiteSizesInUmSq = [160];%In square microns
AnimalMetaData.Probes.SpikeGroupLayoutFilenames = {'NRX_Buzsaki64_5X12';'NRX_Buzsaki64_8X8'};%filenames in /buzcode/GeneralComputation/geometries
AnimalMetaData.Channels.ImpedanceFilenames = {'Ket1_Impedances_5Shank.csv','Ket1_Impedances_8Shank.csv'};%Filenames in basepath

%% Automated after this point
[PerGroupSuperficialToDeep,SpatialXY,NumChansPerProbe,GroupsPerChannel] = bz_ReadProbeMapFiles(AnimalMetaData.Probes.SpikeGroupLayoutFilenames);
AnimalMetaData.Probes.NumGroupsPerProbe = sum(~cellfun(@isempty,PerGroupSuperficialToDeep),2);
AnimalMetaData.Channels.NumChannelsTotal = sum(NumChansPerProbe);
AnimalMetaData.Channels.NumChansPerProbe = NumChansPerProbe;%to do
AnimalMetaData.Channels.ProbeSpikeGroupLayoutSuperficialToDeep = PerGroupSuperficialToDeep;%to do
AnimalMetaData.Channels.ProbeXYLocations = SpatialXY;

%make lookup tables for probe number and anatomy for each channel
po = AnimalMetaData.Probes.PluggingOrder;
lut = [];
lut_ap = [];
glut = [];
glut_ap = [];
pglut = [];
pglut_p = [];
for pidx = 1:AnimalMetaData.Probes.NumberOfProbes
    lut = cat(1,lut,pidx*ones(AnimalMetaData.Channels.NumChansPerProbe(pidx),1));
    lut_ap = cat(1,lut_ap,po(pidx)*ones(AnimalMetaData.Channels.NumChansPerProbe(po(pidx)),1));
    glut = cat(1,glut,GroupsPerChannel{pidx}+length(pglut));
    glut_ap = cat(1,glut_ap,GroupsPerChannel{po(pidx)});
    pglut = cat(1,pglut,[1:AnimalMetaData.Probes.NumGroupsPerProbe(pidx)]'+length(pglut));
    pglut_p = cat(1,pglut_p,pidx*ones(AnimalMetaData.Probes.NumGroupsPerProbe(pidx),1));
end
AnimalMetaData.Probes.GroupsLookupTable = [pglut_p pglut];
AnimalMetaData.Channels.ProbeLookupTable = [[1:AnimalMetaData.Channels.NumChannelsTotal]' lut];
AnimalMetaData.Channels.GroupLookupTable = [[1:AnimalMetaData.Channels.NumChannelsTotal]' glut];
AnimalMetaData.Channels.AnatomyLookupTable = AnimalMetaData.Probes.TargetAnatomy(lut)';
AnimalMetaData.Channels.ProbeLookupTable_AsPlugged = [[1:AnimalMetaData.Channels.NumChannelsTotal]' lut_ap];
AnimalMetaData.Channels.GroupLookupTable_AsPlugged = [[1:AnimalMetaData.Channels.NumChannelsTotal]' glut_ap];
AnimalMetaData.Channels.AnatomyLookupTable_AsPlugged = AnimalMetaData.Probes.TargetAnatomy(lut_ap)';

%Get impedances per channel
AnimalMetaData.Channels.ImpedanceByChannel = cell(AnimalMetaData.Probes.NumberOfProbes);
for pidx = 1:length(AnimalMetaData.Channels.ImpedanceFilenames)
    tf = fullfile(basepath,AnimalMetaData.Channels.ImpedanceFilenames{pidx});
    txt = read_mixed_csv(tf);
    numchans = size(txt,1)-1;
%     channums = 1:numchans;
    for cidx = 1:numchans
        eloc = strfind(txt{cidx+1,5},'e');
        n1 = str2num(txt{cidx+1,5}(1:eloc-1));
        n2 = str2num(txt{cidx+1,5}(eloc+1:end));
        imp(cidx) = n1*10^n2;
    end
%     AnimalMetaData.Probes.ImpedanceByChannel{pidx} = cat(2,channums',imp');
    AnimalMetaData.Channels.ImpedanceByChannel{pidx} = imp';
end

%% someone should combine XY with probe orientation angle and implant 
%coordinates to map each site at day of implant

%% Make XML for animal
aname = AnimalMetaData.AnimalName;
apath = AnimalMetaData.AnimalBaseFolderPath;
pfiles = AnimalMetaData.Probes.SpikeGroupLayoutFilenames;
plugord = AnimalMetaData.Probes.PluggingOrder;
bz_MakeXMLFromProbeMaps(apath,aname,pfiles,plugord);

%% Save
save(fullfile(apath,[aname '.AnimalInfo.mat']),'AnimalMetaData')



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