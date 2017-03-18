function [lfp] = bz_GetLFP(channels,varargin)
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
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
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
	warning('No session defined, so we look for a *lfp file...')
    try
       d = dir('*lfp');
       lfpFile = d.name; % we assume one .lfp file or this should break
       xml = LoadParameters;
       basename = split(lfpFile,'.');
       basename = basename{1};
       path = pwd;
       nChannels = xml.nChannels;
       samplingRate = xml.lfpSampleRate;
    catch
        
    end
else  % backwards compatible with FMAT setcurrentsession
    basename = DATA.session.basename;
    path = DATA.session.path;
    nChannels = DATA.nChannels;
    samplingRate = DATA.rates.lfp;
end

% Default values
intervals = [0 Inf];
select = 'number';

if nargin < 1 | mod(length(varargin),2) ~= 0
  error('Incorrect number of parameters (type ''help <a href="matlab:help GetLFP">GetLFP</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin)
  if ~ischar(varargin{i})
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help GetLFP">GetLFP</a>'' for details).']);
  end
  switch(lower(varargin{i}))
    case {'intervals','restrict'}
      intervals = varargin{i+1};
      if ~isdmatrix(intervals) || size(intervals,2) ~= 2
        error('Incorrect value for property ''intervals'' (type ''help <a href="matlab:help GetLFP">GetLFP</a>'' for details).');
      end
    end
end

filename = [path '/' basename '.lfp'];
if ~exist(filename,'file')
    filename = [path '/' basename '.lfp'];
    if ~exist(filename,'file')
        error(['File ''' filename ''' not found.']);
    end
end
nChannels = nChannels;
if isa(channels,'char') && strcmp(lower(channels),'all')
	channels = (1:nChannels)-1;
end

nIntervals = size(intervals,1);
% returns lfp/bz format
for i = 1:nIntervals
    lfp(i).duration = (intervals(i,2)-intervals(i,1));
    lfp(i).interval = [intervals(i,1) intervals(i,2)];
    
	% Load data and put into struct
	lfp(i).data = LoadBinary(filename,'duration',lfp(i).duration,...
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
        lfp(i).duration = (intervals(i,2)-intervals(i,1));
    end
end

