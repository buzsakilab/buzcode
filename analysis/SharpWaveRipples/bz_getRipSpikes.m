function [spkEventTimes] = bz_getRipSpikes(varargin)
%
%    [SpkEventTimes] = bz_getRipSpikes(varargin)
% Saves spike times of all units inside given events in different ways:
%   1. Absolute and relative time of spikes by unit and by event
%   2. Absolute and relative time of spikes by unit
%   3. Absolute and relative time of spikes by ripple
%
% INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'    full path where session is located (default pwd)
%     'events'      It can be either the following options:
%                   1. Buzcode structure for specific events (ripples, UDStates, ...)
%                      By default it will load ripples (output from bz_DetectSWR).
%                      Specifically, if not provided, it loads this ripple 
%                      structure from 'basepath' (if provided), or from current
%                      folder (if not). The input structure must have the 
%                      following field:
%                        .timestamps: Nx2 matrix with starting and ending times 
%                                     (in sec) of each event.
%                   2. A Nx2 matrix with starting and ending times (in sec)
%                      of each event, just like the .timestamps field.
%                       (N: number of events)
%     'spikes'      buzcode ripple structure (from bz_GetSpikes). 
%                   If not provided, it loads it from 'basepath' (if provided),
%                   or from current folder (if not)
%     'UIDs'        A Mx1 boolean matrix with 1s for units to be considered
%                   and 0s for units to be discarded.(M: number of units)
%     'padding'     additional time after ripple end to still search for spikes. 
%                   (default is 0.05 sec)
%     'saveMat'   	Saves file, logical (default: true) 
%
%    =========================================================================
%
% OUTPUTS
%
% spkEventTimes structure with the followin fields:
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%	  .UnitEventAbs	 MxN cell matrix. In each cell, absolute times of
%                    spikes for that particular unit and event
%	  .UnitEventRel	 MxN cell matrix. In each cell, relative times of
%                    spikes (relative to the starting time of events) for 
%                    that particular unit and event
%	  .UnitAbs		 1xM cell matrix. In each cell, absolute times of
%                    spikes for that particular unit across all events
%	  .UnitRel		 1xM cell matrix. In each cell, relative times of
%                    spikes for that particular unit across all events
%	  .EventAbs		 2xN cell matrix. In the first row, absolute times of
%                    spikes for that particular event across all units. In
%                    the second row, the UID associated to the above spike
%	  .EventRel		 2xN cell matrix. In the first row, relative times of
%                    spikes for that particular event across all units. In
%                    the second row, the UID associated to the above spike
%
%    =========================================================================
%
%    Antonio FR, 2017
%    Converted to buzcode format: Andrea Navas-Olive, 2019

% Parse inputs 
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'events',[], @(x) isnumeric(x) || isstruct(x));
addParameter(p,'spikes',{},@isstruct);
addParameter(p,'UIDs',[],@islogical);
addParameter(p,'padding',0.05,@isnumeric);
addParameter(p,'saveMat', true, @islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
events = p.Results.events;
spikes = p.Results.spikes;
UIDs = p.Results.UIDs;
padding = p.Results.padding;
saveMat = p.Results.saveMat;

% Get session info
basename = bz_BasenameFromBasepath(basepath);
load([basepath filesep basename '.sessionInfo.mat']);

% Default events, UIDs and spikes
if isempty(spikes)
    spikes = load([basepath filesep basename '.spikes.cellinfo.mat']);
    spikes = spikes.spikes;
end
if isempty(UIDs)
    UIDs = ones(size(spikes.UID));
end
if isempty(events)
    events = load([basepath filesep basename '.ripples.events.mat']);
    events = events.ripples;
end
% Starting and ending timestamps
if isnumeric(events)
    timestamps = events;
elseif isstruct(events)
    timestamps = events.timestamps;
else
    warning('Events must be either a Nx2 vector or a bz event structure!');
end

%% Get spikes for each unit and each ripple

% We will save spike times of different units in different ways:
spkEventTimes = {};
% 1. Absolute and relative time of spikes by unit and by ripple
for unit = 1:length(spikes.UID)
    if UIDs(unit)
        for event = 1:length(timestamps)
            % Start and end of ripple
            tini = timestamps(event,1);
            tend = timestamps(event,2) + padding;
            % Spikes of this unit within this ripple interval
            tsUnitEvent = spikes.times{unit};
            tsUnitEvent = tsUnitEvent(tsUnitEvent>=tini & tsUnitEvent<=tend);
            % Absolute time of spikes by unit and by ripple
            spkEventTimes.UnitEventAbs{unit,event} = tsUnitEvent';
            % Relative time of spikes by unit and by ripple to ripple start
            spkEventTimes.UnitEventRel{unit,event} = tsUnitEvent' - tini;
        end
    end
end

% 2. Absolute and relative time of spikes by unit
for unit = 1:length(spikes.UID)
    if UIDs(unit)
        spkEventTimes.UnitAbs{unit} = cell2mat(spkEventTimes.UnitEventAbs(unit,:));
        spkEventTimes.UnitRel{unit} = cell2mat(spkEventTimes.UnitEventRel(unit,:));
    end
end

% 3. Absolute and relative time of spikes by ripple
for event = 1:length(timestamps)
    spkEventTimes.EventAbs{event} = [];
    spkEventTimes.EventRel{event} = [];
    for unit = 1:length(spikes.UID)
        if UIDs(unit)
            spkEventTimes.EventAbs{event} = [ spkEventTimes.EventAbs{event}, ...
                                        [cell2mat(spkEventTimes.UnitEventAbs(unit,event)); ...
                                         cell2mat(spkEventTimes.UnitEventAbs(unit,event))*0+spikes.UID(unit)] ];
            spkEventTimes.EventRel{event} = [ spkEventTimes.EventRel{event}, ...
                                         [cell2mat(spkEventTimes.UnitEventRel(unit,event)); ...
                                          cell2mat(spkEventTimes.UnitEventRel(unit,event))*0+spikes.UID(unit)] ];
        end
    end
    spkEventTimes.EventAbs{event} = sortrows(spkEventTimes.EventAbs{event}')';
    spkEventTimes.EventRel{event} = sortrows(spkEventTimes.EventRel{event}')';
end

% Save
if saveMat
   save(fullfile(basepath, [basename '.spkEventTimes.mat'],'spkEventTimes')); 
end


end
