function parameters = LoadParameters(filename)

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

if nargin < 1 % if we're especially lazy, we assume there is one XML in the current working directory....
   xml = dir('*xml'); 
   filename = xml.name;
end

if ~strcmp(filename(end-3:end),'.xml') % we can now give LoadParameters.m the folder location instead of an actual xml file
    d = dir(fullfile(filename, '*xml'));
    filename = fullfile(filename, d.name);
end

if ~exist(filename),
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
p = convert(t);
% parameters = p;

parameters.session.path = pathname;
parameters.session.name = basename;


if isempty(which('xmltree')),
	error('This function requires the <a href="http://www.artefact.tk/software/matlab/xml/">xmltree</a> toolbox by G. Flandin.');
end


if ~isempty(p.spikeDetection),
	parameters.spikeGroups.nGroups = length(p.spikeDetection.channelGroups.group);
	if parameters.spikeGroups.nGroups == 1,
		parameters.spikeGroups.nSamples = str2num(p.spikeDetection.channelGroups.group.nSamples);
		channels = p.spikeDetection.channelGroups.group.channels.channel;
		if isa(channels,'cell'),
			for channel = 1:length(channels),
				parameters.spikeGroups.groups{1}(channel) = str2num(channels{channel});
			end
		else
			parameters.spikeGroups.groups{1} = str2num(channels);
		end
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
else
	parameters.spikeGroups.nSamples = [];
	parameters.spikeGroups.groups = {};
	parameters.spikeGroups.nGroups = length(p.anatomicalDescription.channelGroups);
end

parameters.nChannels = str2num(p.acquisitionSystem.nChannels);
parameters.nBits = str2num(p.acquisitionSystem.nBits);
parameters.rates.lfp = str2num(p.fieldPotentials.lfpSamplingRate);
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

% for backwards compatibility with loadXml_old.m and variants
% the below code fails with certain XMl files people in the lab use
% so we'll wrap this in a try/catch for now..
try
parameters.Date = p.generalInfo.date;
parameters.VoltageRange = str2num(p.acquisitionSystem.voltageRange);
parameters.Amplification = str2num(p.acquisitionSystem.amplification);
parameters.Offset = str2num(p.acquisitionSystem.offset);
parameters.lfpSampleRate = str2num(p.fieldPotentials.lfpSamplingRate);
catch
   warning('could not load .Date, something may be wrong with your xml...') 
   parameters.lfpSampleRate = parameters.rates.lfp;
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
            for b = 1:length(p.anatomicalDescription.channelGroups.group{a}.channel)
                parameters.AnatGrps(a).Channels(b) = str2num(p.anatomicalDescription.channelGroups.group{a}.channel{b});
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

