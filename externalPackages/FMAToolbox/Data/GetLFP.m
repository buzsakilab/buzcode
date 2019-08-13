function [lfp,indices] = GetLFP(channels,varargin)

%GetLFP - Get local field potentials.
%
%  Load local field potentials from disk (unlike spikes or positions, LFP data
%  is usually too large to keep in memory).
%
%  USAGE
%
%    [lfp,indices] = GetLFP(channels,<options>)
%
%    channels       list of channels to load (use keyword 'all' for all)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'restrict'    list of time intervals to read from the LFP file
%     'intervals'   same as 'restrict' (for backwards compatibility)
%     'select'      select channel by ID ('id', counted from 0 a la NeuroScope)
%                   or by number ('number', counted from 1 a la Matlab)
%                   (default = 'id')
%    =========================================================================
%
%  OUTPUT
%
%    lfp            list of (time,voltage1,...,voltageN) tuples
%    indices        for each tuple, the index of the interval it falls in
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
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


error('This function is now deprecated from /buzcode, try using bz_GetLFP.m')

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Default values
intervals = [0 Inf];
select = 'id';

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help GetLFP">GetLFP</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help GetLFP">GetLFP</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case {'intervals','restrict'},
      intervals = varargin{i+1};
      if ~isdmatrix(intervals) || size(intervals,2) ~= 2,
        error('Incorrect value for property ''intervals'' (type ''help <a href="matlab:help GetLFP">GetLFP</a>'' for details).');
      end
    case 'select',
      select = lower(varargin{i+1});
      if ~isstring_FMAT(select,'id','number'),
        error('Incorrect value for property ''select'' (type ''help <a href="matlab:help GetLFP">GetLFP</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GetLFP">GetLFP</a>'' for details).']);
  end
end

filename = [DATA.session.path '/' DATA.session.basename '.lfp'];
if ~exist(filename,'file'),
    filename = [DATA.session.path '/' DATA.session.basename '.lfp'];
    if ~exist(filename,'file'),
        error(['File ''' filename ''' not found.']);
    end
end
nChannels = DATA.nChannels;
if isa(channels,'char') && strcmp(lower(channels),'all'),
	channels = (1:nChannels)-1;
end

if strcmp(select,'id'),
	channels = channels + 1;
end

nIntervals = size(intervals,1);
lfp = [];
indices = [];
for i = 1:nIntervals,
	duration = (intervals(i,2)-intervals(i,1));
	start = intervals(i,1);
	% Load data
	data = LoadBinary(filename,'duration',duration,'frequency',DATA.rates.lfp,'nchannels',nChannels,'start',start,'channels',channels);
    if isinteger(data)
        error('GetLFP does not work with integer data - use bz_GetLFP instead');
    end
	t = start:(1/DATA.rates.lfp):(start+(length(data)-1)/DATA.rates.lfp);t=t';
	lfp = [lfp ; t data];
	indices = [indices ; i*ones(size(t))];
end

