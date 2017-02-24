%LoadBinary - Load data from a binary file.
%
%  USAGE
%
%    data = LoadBinary(filename,<options>)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'duration'    duration to read (in s) (default = Inf)
%     'frequency'   sampling rate (in Hz) (default = 20kHz)
%     'start'       position to start reading (in s) (default = 0)
%     'nChannels'   number of data channels in the file (default = 1)
%     'channels'    channels to read (default = all)
%     'precision'   sample precision (default = 'int16')
%     'skip'        number of bytes to skip after each value is read
%                   (default = 0)
%    =========================================================================

% Copyright (C) 2004-2006 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function data = LoadBinary(filename,varargin)

% Default values
start = 0;
nChannels = 1;
precision = 'int16';
skip = 0;
duration = Inf;
frequency = 20000;
channels = 1;

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help LoadBinary'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help LoadBinary'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'duration',
      duration = varargin{i+1};
      if ~isa(duration,'numeric') | length(duration) ~= 1 | duration < 0,
        error('Incorrect value for property ''duration'' (type ''help LoadBinary'' for details).');
      end
    case 'frequency',
      frequency = varargin{i+1};
      if ~isa(frequency,'numeric') | length(frequency) ~= 1 | frequency <= 0,
        error('Incorrect value for property ''frequency'' (type ''help LoadBinary'' for details).');
      end
    case 'start',
      start = varargin{i+1};
      if ~isa(start,'numeric') | length(start) ~= 1,
        error('Incorrect value for property ''start'' (type ''help LoadBinary'' for details).');
      end
		if start < 0, start = 0; end
    case 'nchannels',
      nChannels = varargin{i+1};
      if ~IsPositiveInteger(nChannels) | length(nChannels) ~= 1,
        error('Incorrect value for property ''nChannels'' (type ''help LoadBinary'' for details).');
      end
    case 'channels',
      channels = varargin{i+1};
      if ~IsPositiveInteger(channels),
        error('Incorrect value for property ''channels'' (type ''help LoadBinary'' for details).');
      end
    case 'precision',
      precision = varargin{i+1};
      if ~isa(precision,'char'),
        error('Incorrect value for property ''precision'' (type ''help LoadBinary'' for details).');
      end
    case 'skip',
      skip = varargin{i+1};
      if ~IsPositiveInteger(skip) | length(skip) ~= 1,
        error('Incorrect value for property ''skip'' (type ''help LoadBinary'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help LoadBinary'' for details).']);
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

f = fopen(filename,'r');

% Position file index for reading
start = floor(start*frequency)*nChannels*sizeInBytes;
status = fseek(f,start,'bof');
if status ~= 0,
	fclose(f);
	error('Could not start reading (possible reasons include trying to read past the end of the file).');
end

% Determine number of samples when duration is 'inf'

if isinf(duration),
	fileStart = ftell(f);
	status = fseek(f,0,'eof');
	if status ~= 0,
		fclose(f);
		error('Error reading the data file (possible reasons include trying to read past the end of the file).');
	end
	fileStop = ftell(f);
	nSamplesPerChannel = (fileStop-fileStart)/nChannels/sizeInBytes;
	duration = nSamplesPerChannel/frequency;
	frewind(f);
	status = fseek(f,start,'bof');
	if status ~= 0,
		fclose(f);
		error('Could not start reading (possible reasons include trying to read past the end of the file).');
	end
else
	nSamplesPerChannel = floor(frequency*duration);
	if nSamplesPerChannel ~= frequency*duration,
		%disp(['Warning: rounding duration (' num2str(duration) ' -> ' num2str(nSamplesPerChannel/frequency) ')']);
		duration = nSamplesPerChannel/frequency;
	end
end

% For large amounts of data, read chunk by chunk

maxSamplesPerChunk = 100000;
nSamples = nChannels*nSamplesPerChannel;
if nSamples > maxSamplesPerChunk,
	% Determine chunk duration and number of chunks
	nSamplesPerChunk = floor(maxSamplesPerChunk/nChannels)*nChannels;
	durationPerChunk = nSamplesPerChunk/frequency/nChannels;
	nChunks = floor(duration/durationPerChunk);
	% Preallocate memory
	data = zeros(nSamplesPerChannel,length(channels),precision);
	% Read all chunks
	i = 1;
	for j = 1:nChunks,
		d = LoadBinaryChunk(f,'frequency',frequency,'nChannels',nChannels,'channels',channels,'duration',durationPerChunk,'skip',skip);
		[m,n] = size(d);
		if m == 0, break; end
		data(i:i+m-1,:) = d;
		i = i+m;
%  		h=waitbar(j/nChunks);
	end
%  	close(h)
	% If the data size is not a multiple of the chunk size, read the remainder
	remainder = duration - nChunks*durationPerChunk;
	if remainder ~= 0,
		d = LoadBinaryChunk(f,'frequency',frequency,'nChannels',nChannels,'channels',channels,'duration',remainder,'skip',skip);
		[m,n] = size(d);
		if m ~= 0,
			data(i:i+m-1,:) = d;
		end
	end
else
	if skip ~= 0,
		data = fread(f,[nChannels frequency*duration],precision,skip);
	else
		data = fread(f,[nChannels frequency*duration],precision);
	end
	data=data';
	
	if ~isempty(channels),
		data = data(:,channels);
	end
end
fclose(f);
