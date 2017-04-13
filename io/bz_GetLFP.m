function [lfp] = bz_GetLFP(varargin)
%GetLFP - Get local field potentials.
%
%  Load local field potentials from disk. No longer dependent on
%  FMAT/SetCurrentSession.
%
%  USAGE
%
%    [lfp] = GetLFP(channels,<options>)
%
%    channels       list of channels to load (use keyword 'all' for all)
%                   channID is 0-indexing, a la neuroscope
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'fbasename'   base file name to load
%     'restrict'    list of time intervals to read from the LFP file
%     'intervals'   same as 'restrict' (for backwards compatibility)
%    =========================================================================
%
%  OUTPUT
%
%    lfp             struct of lfp data. Can be a single struct or an array
%                    of structs for different intervals.  lfp(1), lfp(2),
%                    etc for intervals(1,:), intervals(2,:), etc
%    .data           [Nt x Nd] matrix of the LFP data
%    .timestamps     [Nt x 1] vector of timestamps to match LFP data
%    .interval       [1 x 2] vector of start/stop times of LFP interval
%    .channels       [Nd X 1] vector of channel ID's
%    .samplingRate   LFP sampling rate [default = 1250]
%    .duration       duration, in seconds, of LFP interval
%
%
%  EXAMPLES
%
%    % channel ID 5 (= # 6), from 0 to 120 seconds
%    lfp = GetLFP(5,'restrict',[0 120]);
%    % same, plus from 240.2 to 265.23 seconds
%    lfp = GetLFP(5,'restrict',[0 120;240.2 265.23]);
%    % multiple channels
%    lfp = GetLFP([1 2 3 4 10 17],'restrict',[0 120]);
%    % channel # 3 (= ID 2), from 0 to 120 seconds
%    lfp = GetLFP(3,'restrict',[0 120],'select','number');

% Copyright (C) 2004-2011 by Michaël Zugaro
% editied by David Tingley, 2017
%
% 'select' option has been removed, it allowed switching between 0 and 1
% indexing.  This should no longer be necessary with .lfp.mat structs
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA)
   p = inputParser;
   addRequired(p,'channels',@isnumeric)
   
   addParameter(p,'fbasename','',@isstr)
   addParameter(p,'restrict',[0 Inf],@isnumeric)
   addParameter(p,'intervals',[0 Inf],@isnumeric)
   
   parse(p,varargin{:})
   fbasename = p.Results.fbasename;
   channels = p.Results.channels;
   % backwards compatible with FMAT
   if sum(p.Results.restrict ~= [0 inf]) > 0
       intervals = p.Results.restrict;
   else
       intervals = p.Results.intervals;
   end
   
   if isempty(fbasename)
       warning('No session defined, so we look for a *lfp file...')
       d = dir('*lfp');
       if length(d) > 1 % we assume one .lfp file or this should break
           error('there is more than one .lfp file in this directory?');
       elseif length(d) == 0
           d = dir('*eeg');
           if isempty(d)
               error('could not find an lfp/eeg file..')
           end
       end
       lfp.Filename = d.name;
       fbasename = strsplit(lfp.Filename,'.');
       fbasename = fbasename{1};
   end
   xml = LoadParameters([fbasename '.xml']);
   path = pwd;
   nChannels = xml.nChannels;
try
    samplingRate = xml.lfpSampleRate;
catch
     samplingRate = xml.rates.lfp; % old ugliness we need to get rid of
end
else  % backwards compatible with FMAT setcurrentsession
   p = inputParser;
   if exist([DATA.session.basename '.lfp'])
       lfp.Filename = [DATA.session.basename '.lfp'];
   elseif exist([DATA.session.basename '.eeg'])
       lfp.Filename = [DATA.session.basename '.eeg']; 
   else
       error('could not find lfp file...')
   end
   
   path = DATA.session.path;
   nChannels = DATA.nChannels;
   samplingRate = DATA.rates.lfp;
   addRequired(p,'channels',@isnumeric)
   addParameter(p,'restrict',[0 Inf],@isnumeric)
   addParameter(p,'intervals',[0 Inf],@isnumeric)
   parse(p,varargin{:})
   channels = p.Results.channels;

   if sum(p.Results.restrict ~= [0 inf]) > 0
       intervals = p.Results.restrict;
   else
       intervals = p.Results.intervals;
   end
end

% we assume 0-indexing like neuroscope, but LoadBinary using 1-indexing to
% load....
channels = channels + 1;


nIntervals = size(intervals,1);
% returns lfp/bz format
for i = 1:nIntervals
    lfp(i).duration = (intervals(i,2)-intervals(i,1));
    lfp(i).interval = [intervals(i,1) intervals(i,2)];
    
	% Load data and put into struct
	lfp(i).data = LoadBinary(lfp.Filename,'duration',lfp(i).duration,...
                  'frequency',samplingRate,'nchannels',nChannels,...
                  'start',lfp(i).interval(1),'channels',channels);
	lfp(i).timestamps = [lfp(i).interval(1):(1/samplingRate):...
                        (lfp(i).interval(1)+(length(lfp(i).data)-1)/...
                        samplingRate)]';
    lfp(i).channels = channels;
    lfp(i).samplingRate = samplingRate;
    % check if duration is inf, and reset to actual duration...
    if lfp(i).interval(2) == inf
        lfp(i).interval(2) = length(lfp(i).timestamps)/lfp(i).samplingRate;
        lfp(i).duration = (lfp(i).interval(i,2)-lfp(i).interval(i,1));
    end
end

