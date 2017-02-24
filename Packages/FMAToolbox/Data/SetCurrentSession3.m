function SetCurrentSession3(filename)

%SetCurrentSession3 - Load all data for a given recording session.
%
% This is a modified version of SetCurrentSession() that:
%
% 1) defaults to Spotter's default video parameters when no video files are present.
% 2) does NOT load in spike files (which saves time if all you want is the
% other info
%
% - Luke Sjulson, 2014-08-13
%
% Set current session files and read data from disk. Calling SetCurrentSession
% without parameters will display a file selection dialog.
%
%  USAGE
%
%    SetCurrentSession(filename)
%
%    filename            optional parameter file name; use 'same' to force reload
%
%  NOTE
%
%    If no parameter file name is specified, an interactive file selection
%    dialog is displayed.

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
separator = filesep;

% Initialization
if isempty(DATA) || ~isfield(DATA,'session') || ~isfield(DATA.session,'path') || ~isfield(DATA.session,'basename'),
	format long g;
	DATA.session.basename = '';
	DATA.session.path = '';
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
	Settings;
end

if nargin == 0 || (strcmp(filename,'same') && isempty(DATA.session.basename)),
	% Interactive mode
	[filename,path] = uigetfile('*.xml','Please select a parameter file for this session');
	if filename == 0,return; end
	filename = [path filename];
end

if strcmp(filename,'same'),
	% Force reload
	path = DATA.session.path;
	basename = DATA.session.basename;
else
	% Parse file name
	[path,basename] = fileparts(filename);
	if isempty(path)|path(1)~=separator, path = [pwd separator path]; end
	path = strrep(path,[separator '.' separator],separator);
end

disp(['Loading session files for ' basename]);

% File already loaded?
if strcmp(basename,DATA.session.basename) & strcmp(path,DATA.session.path) & ~strcmp(filename,'same'),
	disp(['... session files already loaded, skipping - type SetCurrentSession(''same'') to force reload']);
	disp('Done');
	return
end

% Parameter file
%DATA = LoadParameters([path separator basename '.xml']);
DATA = LoadParameters2([basename '.xml']);
disp(['... loaded parameter file ''' basename '.xml''']);

% Event file(s)
DATA.events.time = [];
DATA.events.description = {};
eventFiles = dir([path separator basename '.*.evt']);
if ~isempty(eventFiles),
	for i = 1:length(eventFiles),
		events = LoadEvents([path separator eventFiles(i).name]);
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
if exist([path separator basename '.pos']),
	DATA.positions = LoadPositions([path separator basename '.pos'],DATA.rates.video);
	disp(['... loaded position file ''' basename '.pos''']);
elseif exist([path separator basename '.whl']),
	DATA.positions = LoadPositions([path separator basename '.whl'],DATA.rates.video);
	disp(['... loaded position file ''' basename '.whl''']);
elseif exist([path separator basename '.whl']),
	DATA.positions = LoadPositions([path separator basename '.mqa'],DATA.rates.video);
	disp(['... loaded position file ''' basename '.mqa''']);
else
	disp('... (no position file found)');
end

% % Spike files
% DATA.spikes = [];
% for i = 1:DATA.spikeGroups.nGroups,
% 	filename = [path separator basename '.' int2str(i) '.clu'];
% 	if exist(filename,'file'),
% 		try
% 			DATA.spikes = [DATA.spikes;LoadSpikeTimes(filename,DATA.rates.wideband)];
% 			disp(['... loaded spike files ''' basename '.' int2str(i) '.clu''']);
% 		catch
% 			disp(['... (could not load spike files ''' basename '.' int2str(i) '.clu'')']);
% 		end
% 	else
% 		filename = [path separator basename '.clu.' int2str(i)];
% 		if exist(filename,'file'),
% 			try
% 				DATA.spikes = [DATA.spikes;LoadSpikeTimes(filename,DATA.rates.wideband)];
% 				disp(['... loaded spike files ''' basename '.clu.' int2str(i) '''']);
% 			catch
% 				disp(['... (could not load spike files ''' basename '.clu.' int2str(i) ''')']);
% 			end
% 		end
% 	end
% end
% if isempty(DATA.spikes),
% 	disp('... (no spike files found)');
% end

% This is updated only once the files have been properly loaded
DATA.session.basename = basename;
DATA.session.path = path;

disp('Done');

