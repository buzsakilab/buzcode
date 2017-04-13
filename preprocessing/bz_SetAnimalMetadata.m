function AnimalMetadata = bz_SetAnimalMetadata(basepath,basename)
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
%    baseName.AnimalMetadata.mat written to basepath
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

%% Initial variable parsing
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

%% Value setting - humans must do this.  Will be asked to edit a _NoteText.m file
notesname = [basename,'_AnimalNotesText.m'];
if ~exist(fullfile(basepath,notesname),'file')
    w = which('bz_AnimalNotesTemplate.m');% copy an example header here to edit
    copyfile(w,notesname);
end

edit(notesname)
prompt = 'Push any key in this window when done editing the NoteText file ';
str = input(prompt,'s');
run(notesname);%save _AnimalNotes.mat to disk

load(fullfile(basepath,[basename '_AnimalNotes.mat']))%load AnimalNotes
AnimalMetadata = AnimalNotes;
clear AnimalNotes

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

%% Make an initial .xml file
pfiles = AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames;
plugord = AnimalMetadata.ExtracellEphys.Probes.PluggingOrder;
dfs = AnimalMetadata.ExtracellEphys.Parameters;

xmlfilename = fullfile(basepath,[basename,'.xml']);
if ~exist(xmlfilename,'file')
    bz_MakeXMLFromProbeMaps(basepath,basename,pfiles,plugord,dfs);
else
    display([basename,'.xml is already in existence, please rename or ',...
        'delete if you would like to overwrite an existing .xml'])
end


%% Save
save(fullfile(basepath,[basename '_AnimalMetadata.mat']),'AnimalMetadata')



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