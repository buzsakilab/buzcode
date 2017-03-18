function parameters = LoadParameters(filename)

%LoadParameters - Load parameters from an XML file.
%
%  USAGE
%
%    parameters = LoadParameters(filename)
%
%    filename            parameter file name

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% updated for compatibility by David Tingley 02/2017

if nargin < 1 % if we're especially lazy, we assume there is one XML in the current working directory....
   xml = dir('*xml'); 
   filename = xml.name;
end

if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
[pathname,basename] = fileparts(filename);

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
	parameters.spikeGroups.nGroups = 0;
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
parameters.nElecGps = length(p.spikeDetection.channelGroups.group);
parameters.ElecGp = p.spikeDetection.channelGroups.group;
parameters.HiPassFreq = 500; % default hi-pass for klusta-3.0 w/ intan data

% for backwards compatibility with loadXml_old.m and variants
parameters.Date = p.generalInfo.date;
parameters.VoltageRange = str2num(p.acquisitionSystem.voltageRange);
parameters.Amplification = str2num(p.acquisitionSystem.amplification);
parameters.Offset = str2num(p.acquisitionSystem.offset);
parameters.lfpSampleRate = str2num(p.fieldPotentials.lfpSamplingRate);

% fixing AnatGrps and SpkGrps
for a = 1:length(p.anatomicalDescription.channelGroups.group)
    for b = 1:length(p.anatomicalDescription.channelGroups.group{a}.channel)
        parameters.AnatGrps(a).Channels(b) = str2num(p.anatomicalDescription.channelGroups.group{a}.channel{b});
    end        
end
for a = 1:parameters.spikeGroups.nGroups
    parameters.SpkGrps(a).Channels = parameters.spikeGroups.groups{a};
    parameters.SpkGrps(a).nSamples =  str2num(p.spikeDetection.channelGroups.group{a}.nSamples);
    parameters.SpkGrps(a).PeakSample = str2num(p.spikeDetection.channelGroups.group{a}.peakSampleIndex); 
    parameters.SpkGrps(a).nFeatures =  str2num(p.spikeDetection.channelGroups.group{a}.nFeatures); 
end

