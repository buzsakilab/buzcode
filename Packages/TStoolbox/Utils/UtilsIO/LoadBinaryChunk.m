%LoadBinaryChunk - Load data chunck from an open binary file.
%
%  USAGE
%
%    data = LoadBinaryChunk(fid,<options>)
%
%    fid            file id (obtained via fopen)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'duration'    duration to read (in s) (default = 1s)
%     'frequency'   sampling rate (in Hz) (default = 20kHz)
%     'start'       position to start reading (in s) (default = from current
%                   index, allowing to read a file chunck by chunck)
%     'nChannels'   number of data channels in the file (default = 1)
%     'precision'   sample precision (default = 'int16')
%    =========================================================================

% Copyright (C) 2004-2006 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function data = LoadBinaryChunk(fid,varargin)

% Default values
start = 0;
fromCurrentIndex = true;
nChannels = 1;
precision = 'int16';
duration = 1;
frequency = 20000;
channels = [];

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help LoadBinaryChunk'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help LoadBinaryChunk'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'duration',
      duration = varargin{i+1};
      if ~isa(duration,'numeric') | length(duration) ~= 1 | duration < 0,
        error('Incorrect value for property ''duration'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'frequency',
      frequency = varargin{i+1};
      if ~isa(frequency,'numeric') | length(frequency) ~= 1 | frequency <= 0,
        error('Incorrect value for property ''frequency'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'start',
      start = varargin{i+1};
      fromCurrentIndex = false;
      if ~isa(start,'numeric') | length(start) ~= 1,
        error('Incorrect value for property ''start'' (type ''help LoadBinaryChunk'' for details).');
      end
		if start < 0, start = 0; end
    case 'nchannels',
      nChannels = varargin{i+1};
      if ~IsPositiveInteger(nChannels) | length(nChannels) ~= 1,
        error('Incorrect value for property ''nChannels'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'channels',
      channels = varargin{i+1};
      if ~IsPositiveInteger(channels),
        error('Incorrect value for property ''channels'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'precision',
      precision = varargin{i+1};
      if ~isa(precision,'char'),
        error('Incorrect value for property ''precision'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'skip',
      skip = varargin{i+1};
      if ~IsPositiveInteger(skip) | length(skip) ~= 1,
        error('Incorrect value for property ''skip'' (type ''help LoadBinaryChunk'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help LoadBinaryChunk'' for details).']);
  end
end

sizeInBytes = 0;
switch precision,
	case {'uchar','unsigned char','schar','signed char','int8','integer*1','uint8','integer*1'},
		sizeInBytes = 1;
	case {'int16','integer*2','uint16','integer*2'},
		sizeInBytes = 2;
	case {'int32','integer*4','uint32','integer*4','single','real*4','float32','real*4'},
		sizeInBytes = 4;
	case {'int64','integer*8','uint64','integer*8','double','real*8','float64','real*8'},
		sizeInBytes = 8;
end

% Position file index for reading
if ~fromCurrentIndex,
	start = floor(start*frequency)*nChannels*sizeInBytes;
	status = fseek(fid,start,'bof');
	if status ~= 0,
		error('Could not start reading (possible reasons include trying to read a closed file or past the end of the file).');
	end
end

% Read data chunck
if skip ~= 0,
	data = fread(fid,[nChannels frequency*duration],precision,skip);
else
	data = fread(fid,[nChannels frequency*duration],precision);
end;
data=data';

% Keep only required channels
if ~isempty(channels) & ~isempty(data),
	data = data(:,channels);
end

