function [parameters] = LoadParameters(filename)

%LoadParameters - Load parameters from an XML file.
%
%  USAGE
%
%    parameters = LoadParameters(filename)
%
%    filename            parameter file name
%       -or-
%    basePath            where the parameter file baseName.xml lives

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% updated for compatibility by David Tingley 02/2017
% updated for compatibility by Rachel Swanson 05/28/2017
% updated for compatibility by Daniel Levenstein 09/2017

if nargin < 1 % if we're especially lazy, we assume there is one XML in the current working directory....
%    xml = dir('*xml'); 
%    filename = xml.name;
    filename = pwd; %if no input filename is basepath and you're looking for basename.xml
end

if ~strcmp(filename(end-3:end),'.xml') % we can now give LoadParameters.m the folder location instead of an actual xml file
    d = dir(fullfile(filename, '*xml')); 
    
    if length(d)>1 %if multiple .xmls, pick the one that matches baseName
        baseName = bz_BasenameFromBasepath(filename);
        correctxmlfilename = cellfun(@(X) strcmp(X,[baseName,'.xml']),{d.name});
        d = d(correctxmlfilename);
        display(['Multiple .xml files in this folder, trying ',baseName,'.xml'])
    end
    if isempty(d) %if no .xmls - you have a problem
        error('LoadParameters:noXmls',['No .xml in ',filename])
    end
    
    filename = fullfile(filename, d.name);
end

if ~exist(filename)
	error(['File ''' filename ''' not found.']);
end
[pathname,basename,extension] = fileparts(filename);
if ~strcmp(extension,'.xml')
    basename = [basename extension];
end

if isempty(pathname)
    pathname = pwd;
end

t = xmltree(filename);

paths = which('convert','-ALL');
for i = 1:length(paths) % we need to find the right version of convert.m (2018a bugfix)
    if ~isempty(strfind(paths{i},'buzcode'))
       addpath(paths{i},'-begin') % bump buzcode/XML/convert to the top of the search path
    end
end
p = convert(t);
% parameters = p;

parameters.session.path = pathname;
parameters.session.name = basename;


if isempty(which('xmltree')),
	error('This function requires the <a href="http://www.artefact.tk/software/matlab/xml/">xmltree</a> toolbox by G. Flandin.');
end


if ~isempty(p.spikeDetection),
    try
	parameters.spikeGroups.nGroups = length(p.spikeDetection.channelGroups.group);
	if parameters.spikeGroups.nGroups == 1,   %if there's a single spike group		
        if isfield(p.spikeDetection.channelGroups.group,'nSamples')
            parameters.spikeGroups.nSamples = str2num(p.spikeDetection.channelGroups.group.nSamples);
        end
		channels = p.spikeDetection.channelGroups.group.channels.channel;
		if isa(channels,'cell'),
			for channel = 1:length(channels),
				parameters.spikeGroups.groups{1}(channel) = str2num(channels{channel});
			end
		else
			parameters.spikeGroups.groups{1} = str2num(channels);
        end
        
        %%%%%%%bug fix for xml tree mishandling single spike group
         temp = p.spikeDetection.channelGroups.group;
         p.spikeDetection.channelGroups = rmfield(p.spikeDetection.channelGroups,'group');
         p.spikeDetection.channelGroups.group{1} = temp;  
	else
		for group = 1:parameters.spikeGroups.nGroups,
            if isfield(p.spikeDetection.channelGroups.group{group},'nSamples')
    			parameters.spikeGroups.nSamples(group) = str2num(p.spikeDetection.channelGroups.group{group}.nSamples);
            end
			channels = p.spikeDetection.channelGroups.group{group}.channels.channel;
			if isa(channels,'cell'),
				for channel = 1:length(channels),
					parameters.spikeGroups.groups{group}(channel) = str2num(channels{channel});
				end
			else
				parameters.spikeGroups.groups{group} = str2num(channels);
			end
		end
    end
    catch
        warning('something went wrong loading spikeGroups from XML')
        parameters.spikeGroups.nSamples = [];
        parameters.spikeGroups.groups = {};
        parameters.spikeGroups.nGroups = 0;
    end
else
	parameters.spikeGroups.nSamples = [];
	parameters.spikeGroups.groups = {};
	parameters.spikeGroups.nGroups = 0;
end

parameters.nChannels = str2num(p.acquisitionSystem.nChannels);
%channel list assumes 0-indexing a la neuroscope, and that all channels are used.
%In the future, should think about this...
parameters.channels = [0:parameters.nChannels-1]; 
parameters.nBits = str2num(p.acquisitionSystem.nBits);

% Use manually-entered LFP sample rate if user entered it in Neuroscope
manuallfpsamprate = 0;
if isfield(p,'files')
    if isfield(p.files,'file')
        if isfield(p.files.file,'extension')
            if strcmp(p.files.file.extension,'lfp')
                manuallfpsamprate = 1;
            end
        end
    end
end
if manuallfpsamprate
    parameters.rates.lfp = str2num(p.files.file.samplingRate);
else
    parameters.rates.lfp = str2num(p.fieldPotentials.lfpSamplingRate);
end            
parameters.rates.wideband = str2num(p.acquisitionSystem.samplingRate);
try
	parameters.rates.video = str2num(p.video.samplingRate);
% 	parameters.maxX = str2num(p.video.width);   % Deprecated with new tracking systems
% 	parameters.maxY = str2num(p.video.height);  % Deprecated with new tracking systems
catch
	parameters.rates.video = 0;
% 	parameters.maxX = 0;  % Deprecated with new tracking systems
% 	parameters.maxY = 0;  % Deprecated with new tracking systems 
% 	disp('... warning: missing video parameters (set to zero)');
end

% for backwards compatibility with LoadPar.m

parameters.FileName = parameters.session.name;  % killing me slowly with redundancy 
parameters.SampleTime = (1/str2num(p.acquisitionSystem.samplingRate)) * 1e+6; % duration of a sample in microseconds
parameters.nElecGps = length(p.anatomicalDescription.channelGroups.group);
parameters.ElecGp = p.anatomicalDescription.channelGroups.group;
parameters.HiPassFreq = 500; % default hi-pass for klusta-3.0 w/ intan data
parameters.lfpSampleRate = parameters.rates.lfp;

% for backwards compatibility with loadXml_old.m and variants
% the below code fails with certain XMl files people in the lab use
% so we'll wrap this in a try/catch for now..
try
    parameters.VoltageRange = str2num(p.acquisitionSystem.voltageRange);
    parameters.Amplification = str2num(p.acquisitionSystem.amplification);
    parameters.Offset = str2num(p.acquisitionSystem.offset);
catch
     warning('could not load .Amplification, something may be wrong with your xml...') 
end
try
parameters.Date = p.generalInfo.date;
% 
% catch
%    warning('could not load .Date, something may be wrong with your xml...') 
end
% fixing AnatGrps and SpkGrps
% the below code fails with certain XMl files people in the lab use
% so we'll wrap this in a try/catch for now..
try
    for a = 1:length(p.anatomicalDescription.channelGroups.group)
        if isstruct(p.anatomicalDescription.channelGroups.group)
            for b = 1:length(p.anatomicalDescription.channelGroups.group(a).channel)
                parameters.AnatGrps(a).Channels(b) = str2num(p.anatomicalDescription.channelGroups.group(a).channel{b});
            end 
        elseif iscell(p.anatomicalDescription.channelGroups.group)
            if iscell(p.anatomicalDescription.channelGroups.group{a}.channel)
                for b = 1:length(p.anatomicalDescription.channelGroups.group{a}.channel)
                    parameters.AnatGrps(a).Channels(b) = str2num(p.anatomicalDescription.channelGroups.group{a}.channel{b});
                end 
            elseif isvector(p.anatomicalDescription.channelGroups.group{a}.channel)
                parameters.AnatGrps(a).Channels = str2num(p.anatomicalDescription.channelGroups.group{a}.channel);
            else
                warning('Anatomy Groups seems to have an issue, eh?..') 
            end

        end
    end
    for a = 1:parameters.spikeGroups.nGroups
        if ~isempty(parameters.spikeGroups.nSamples)
        parameters.SpkGrps(a).Channels = parameters.spikeGroups.groups{a};
        parameters.SpkGrps(a).nSamples =  str2num(p.spikeDetection.channelGroups.group{a}.nSamples);
        parameters.SpkGrps(a).PeakSample = str2num(p.spikeDetection.channelGroups.group{a}.peakSampleIndex); 
        parameters.SpkGrps(a).nFeatures =  str2num(p.spikeDetection.channelGroups.group{a}.nFeatures); 
        else
            for b = 1:length(p.anatomicalDescription.channelGroups.group.channel)
                parameters.SpkGrps(a).Channels(b) = str2num(p.anatomicalDescription.channelGroups.group.channel{b});
                parameters.SpkGrps(a).nSamples =  [];
                parameters.SpkGrps(a).PeakSample = [];
                parameters.SpkGrps(a).nFeatures =  [];  
            end
        end
    end
catch
  warning('could not load .SpkGrps and .AnatGrps, something may be missing from your XML file..') 
end

%% Unit Info - will generally just get passed through in bz_GetSpikes
try
    numunits = length(p.units.unit);
    for uu = 1:numunits
        parameters.Units(uu).spikegroup = str2num(p.units.unit{uu}.group);
        parameters.Units(uu).cluster = str2num(p.units.unit{uu}.cluster);
        parameters.Units(uu).structure = p.units.unit{uu}.structure;
        parameters.Units(uu).type = p.units.unit{uu}.type;
        parameters.Units(uu).isolationDistance = str2num(p.units.unit{uu}.isolationDistance);
        parameters.Units(uu).quality = p.units.unit{uu}.quality;
        parameters.Units(uu).notes = p.units.unit{uu}.notes;
    end
catch
%if no units in the .xml..... well probably not using them anyway, eh?
end
%% For added plugins (such as badchannels, regions)
try %some xml may not have p.programs.program.... if so, ignore all of this
    plugins = p.programs.program;
    pluginnames = cellfun(@(X) X.name,plugins,'uniformoutput',false);
    %Run through each plugin and check if it matches something we know what
    %to do with, feel free to add more things here for your own purposes
    for pp = 1:length(pluginnames)
        
        if strcmp(pluginnames{pp},'badchannels')
            %Badchannels should be a plugin in the xml, with a single
            %parameter "badchannels" and a list of bad channels separated
            %by a space. This is temporary while we get metadata ironed out
            %-DL
            assert(strcmp(plugins{pp}.parameters.parameter.name,'badchannels'),...
                'There is a plugin ''badchannels'', but the parameter name is not ''badchannels''')
            parameters.badchannels = str2num(plugins{pp}.parameters.parameter.value);
        end
        
        if strcmp(pluginnames{pp},'regions')
            %Regions should be a plugin in the xml, with a group per region
            %with the region name and a list of channels separated
            %by a space or comma. This is temporary while we get metadata ironed out
            %-DL
            parameters.region = repmat({''},1,sessionInfo.nChannels);
            numregions = length(plugins{pp}.parameters.parameter);
            for rr = 1:numregions
                regionname = plugins{pp}.parameters.parameter{rr}.name;
                if strcmp(regionname,'regionname') %"default" is 'regionname'
                    continue
                end
                regionchans = str2num(plugins{pp}.parameters.parameter{rr}.value);
                parameters.region(ismember(parameters.channels,regionchans)) = {regionname};
            end
        end
    end
end

