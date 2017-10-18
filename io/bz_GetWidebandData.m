function [data,indices] = bz_GetWidebandData(channels,varargin)

% bz_GetWidebandData - Get local field potentials.
%
%  Load wide band data from disk (unlike spikes or positions, raw data is way
%  too large to keep in memory).
%
%  USAGE
%
%    [data,indices] = bz_GetWidebandData(channels,<options>)
%
%    channels       optional list of channels to load (default = all)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'restrict'    list of time intervals to read from the data file
%     'intervals'   same as 'restrict' (for backwards compatibility)
%     'select'      select channel by ID ('id', counted from 0 a la NeuroScope)
%                   or by number ('number', counted from 1 a la Matlab)
%                   (default = 'id')
%    =========================================================================
%
%  OUTPUT
%
%    data           list of (time,voltage1,...,voltageN) tuples
%    indices        for each tuple, the index of the interval it falls in
%
%  EXAMPLES
%
%    % channel ID 5 (= # 6), from 0 to 120 seconds
%    data = bz_GetWidebandData(5,'intervals',[0 120]);
%    % same, plus from 240.2 to 265.23 seconds
%    data = bz_GetWidebandData(5,'intervals',[0 120;240.2 265.23]);
%    % multiple channels
%    data = bz_GetWidebandData([1 2 3 4 10 17],'intervals',[0 120]);
%    % channel # 3 (= ID 2), from 0 to 120 seconds
%    data = bz_GetWidebandData(3,'intervals',[0 120],'select','number');

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Default values
intervals = [0 Inf];
select = 'id';

% Optional parameter
if ischar(channels),
	varargin = {channels,varargin{:}};
	channels = []; % all
end

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help bz_GetWidebandData">bz_GetWidebandData</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help bz_GetWidebandData">bz_GetWidebandData</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case {'intervals','restrict'},
      intervals = varargin{i+1};
      if ~isdmatrix(intervals) || size(intervals,2) ~= 2,
        error('Incorrect value for property ''intervals'' (type ''help <a href="matlab:help bz_GetWidebandData">bz_GetWidebandData</a>'' for details).');
      end
    case 'select',
      select = lower(varargin{i+1});
      if ~isstring_FMAT(select,'id','number'),
        error('Incorrect value for property ''select'' (type ''help <a href="matlab:help bz_GetWidebandData">bz_GetWidebandData</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help bz_GetWidebandData">bz_GetWidebandData</a>'' for details).']);
  end
end

filename = [DATA.session.path '/' DATA.session.basename '.dat'];
nChannels = DATA.nChannels;
if isempty(channels),
	channels = 1:nChannels;
elseif strcmp(select,'id'),
	channels = channels + 1;
end

nIntervals = size(intervals,1);
data = [];
indices = [];
for i = 1:nIntervals,
	duration = (intervals(i,2)-intervals(i,1));
	start = intervals(i,1);
	% Load data
	d = bz_LoadBinary(filename,'duration',duration,'frequency',DATA.rates.wideband,'nchannels',nChannels,'start',start,'channels',channels);
	% The following two lines compensate for annoying numerical precision errors in Matlab, whereby the number of samples
	% read with LoadBinary is not always exactly that expected, depending on how accurately 'duration' is coded internally
	n = size(d,1);
	t = linspace(start,start+n/DATA.rates.wideband,n)';
	data = [data ; t d];
	indices = [indices ; i*ones(size(t))];
end

