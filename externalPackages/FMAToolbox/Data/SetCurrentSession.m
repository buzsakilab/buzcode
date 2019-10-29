function SetCurrentSession(varargin)

%SetCurrentSession - Load all data for a given recording session.
%
% Set current session files and read data from disk. Calling SetCurrentSession
% without parameters will display a file selection dialog.
%
%  USAGE
%
%    SetCurrentSession(varargin)
%
%    basepath       optional path to recording you'd like to load
%    basename       optional recording name to load
%    xmlFile        optional xml filename to load
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'spikes'      logical to load spikes (default = true)
%    =========================================================================
%
%  EXAMPLES
%   
%   SetCurrentSession('basepath',pwd)
%   SetCurrentSession('basepath',pwd,'basename','20170313_0um_0um_170313_122027','xmlFile','20170313_0um_0um_170313_122027.xml')
%
%   Can also be used without basepath but you MUST BE IN THE RECORDING DIRECTORY
%
%   SetCurrentSession
%   SetCurrentSession('basename','20170313_0um_0um_170313_122027')  
%   SetCurrentSession('xmlFile','20170313_0um_0um_170313_122027.xml')
%   SetCurrentSession('basename','20170313_0um_0um_170313_122027','xmlFile','20170313_0um_0um_170313_122027.xml') 
%
%
%  NOTE
%
%    If no xmlFile is givien, it is first assumed to be the same naming as 'basename'
%    If if that xmlFile can't be found, an interactive file selection
%    dialog is displayed.
%    If an xmlFile is given that does not exist, an error will be thrown

% Copyright (C) 2004-2014 by Michaël Zugaro, 2014 by Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


p = inputParser;
addParameter(p,'basepath','',@isstr)
addParameter(p,'basename','',@isstr)
addParameter(p,'xmlFile','',@isstr)
addParameter(p,'spikes',true,@islogical)
addParameter(p,'reload',false,@islogical)

% addParameter(p,'position',true,@islogical)  to be added some day
% addParameter(p,'events',true,@islogical)

parse(p,varargin{:})

basename = p.Results.basename;
basepath = p.Results.basepath;
xmlFile = p.Results.xmlFile;
spikes = p.Results.spikes;
reload = p.Results.reload;

% did you give one of the correct input sets? let's find out...

% if you gave nothing, you better be in the right directory....
if isempty(basepath) & isempty(basename) & isempty(xmlFile)
   basepath = pwd;
   temp = strsplit(basepath,filesep);
   basename = temp{length(temp)};
   xmlFile = [basename '.xml'];
   if ~exist([basepath filesep xmlFile])
      [xmlFile,basepath] = uigetfile('*.xml','Please select a parameter file for this session'); 
      basename = xmlFile(1:end-3);
   end
end

% if you gave only basepath, there better be an xml there...
if ~isempty(basepath) & isempty(basename) & isempty(xmlFile)
   temp = strsplit(basepath,filesep);
   basename = temp{length(temp)};
   xmlFile = [basename '.xml'];
   if ~exist([basepath filesep xmlFile])
      error('could not find an xml...') 
   end
end

% if you gave basepath and xmlFile, let's check if that xml exists...
if ~isempty(basepath) & isempty(basename) & ~isempty(xmlFile)
   temp = strsplit(basepath,filesep);
   basename = temp{length(temp)};
   if ~exist([basepath filesep xmlFile])
      error('could not find an xml...') 
   end
end

% if you gave only basename, you better be in the right directory....
if isempty(basepath) & ~isempty(basename) & isempty(xmlFile)
   basepath = pwd;
   xmlFile = [basename '.xml'];
   if ~exist([basepath filesep xmlFile])
      error('could not find an xml...') 
   end
end

% if you gave only xmlFile, you better be in the right directory....
if isempty(basepath) & isempty(basename) & ~isempty(xmlFile)
   basepath = pwd;
   temp = strsplit(basepath,filesep);
   basename = temp{length(temp)};
   if ~exist([basepath filesep xmlFile])
      error('could not find an xml...') 
   end
end

global DATA;

% Initialization
if isempty(DATA) || ~isfield(DATA,'session') || ~isfield(DATA.session,'basepath') || ~isfield(DATA.session,'basename')
	format long g;
	DATA.session.basename = [];
	DATA.session.basepath = [];
	DATA.spikeGroups.nGroups = 0;
	DATA.spikeGroups.nSamples = [];
	DATA.spikeGroups.peakSamples = [];
	DATA.spikeGroups.groups = {};
	DATA.nChannels = [];
	DATA.nBits = [];
	DATA.rates.lfp = [];
	DATA.rates.wideband = [];
	DATA.rates.video = [];
	DATA.maxX = [];
	DATA.maxY = [];
	DATA.events.time = [];
	DATA.events.description = {};
	DATA.positions = [];
	DATA.spikes = [];
	% Default settings
%  	Settings;
end

disp(['Loading session files for ' basename]);

% File already loaded?
if strcmp(basename,DATA.session.basename) && strcmp(basepath,DATA.session.basepath) && reload == 0
	disp(['... session files already loaded, skipping - type SetCurrentSession(''same'') to force reload']);
	disp('Done');
	return
end

% Parameter file
DATA = LoadParameters([basepath filesep xmlFile]);
disp(['... loaded parameter file ''' xmlFile '']);

% Event file(s)
DATA.events.time = [];
DATA.events.description = {};
eventFiles = dir([basepath filesep basename '.*.evt']);
if ~isempty(eventFiles),
	for i = 1:length(eventFiles),
		events = LoadEvents([basepath filesep eventFiles(i).name]);
		if isempty(events.time), continue; end
		DATA.events.time = [DATA.events.time ; events.time];
		DATA.events.description = {DATA.events.description{:} events.description{:}}';
		disp(['... loaded event file ''' eventFiles(i).name '''']);
	end
	[DATA.events.time,i] = sortrows(DATA.events.time);
	DATA.events.description = {DATA.events.description{i}}';
else
	disp('... (no event file found)');
end

% Position file
DATA.positions = [];
if exist([basepath filesep basename '.pos']),
	DATA.positions = LoadPositions([basepath filesep basename '.pos'],DATA.rates.video);
	disp(['... loaded position file ''' basename '.pos''']);
elseif exist([basepath filesep basename '.whl']),
	DATA.positions = LoadPositions([basepath filesep basename '.whl'],DATA.rates.video);
	disp(['... loaded position file ''' basename '.whl''']);
elseif exist([basepath filesep basename '.whl']),
	DATA.positions = LoadPositions([basepath filesep basename '.mqa'],DATA.rates.video);
	disp(['... loaded position file ''' basename '.mqa''']);
else
	disp('... (no position file found)');
end

% Spike files
if spikes == 1
	DATA.spikes = [];
	for i = 1:DATA.spikeGroups.nGroups,
		fullpath = [basepath filesep basename '.' int2str(i) '.clu'];
		if exist(fullpath,'file'),
			try
				DATA.spikes = [DATA.spikes;LoadSpikeTimes(fullpath,DATA.rates.wideband)];
				disp(['... loaded spike files ''' fullpath '.' int2str(i) '.clu''']);
			catch
				disp(['... (could not load spike files ''' fullpath '.' int2str(i) '.clu'')']);
			end
		else
			fullpath = [basepath filesep basename '.clu.' int2str(i)];
			if exist(fullpath,'file'),
				try
					DATA.spikes = [DATA.spikes;LoadSpikeTimes(fullpath,DATA.rates.wideband)];
					disp(['... loaded spike files ''' fullpath '.clu.' int2str(i) '''']);
				catch
					disp(['... (could not load spike files ''' fullpath '.clu.' int2str(i) ''')']);
				end
			end
		end
	end
	if isempty(DATA.spikes),
		disp('... (no spike files found)');
	end
else
	disp('... (skipping spike files)');
end

% This is updated only once the files have been properly loaded
DATA.session.basename = basename;
DATA.session.basepath = basepath;

disp('Done');
