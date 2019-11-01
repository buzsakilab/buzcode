function ChangeBinaryGain(filename,varargin)

%ChangeBinaryGain - Change gain in a multiplexed binary file.
%
%  USAGE
%
%    ChangeBinaryGain(filename,<options>)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'output'      output file name (default = append '-gains')
%     'nChannels'   number of data channels in the file (default = 1)
%     'gains'       list of gains for each channel, or one value for all
%     'precision'   sample precision (default = 'int16')
%     'skip'        number of bytes to skip after each value is read
%                   (default = 0)
%    =========================================================================

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
nChannels = 1;
precision = 'int16';
warning(['this function assumes int16 precision, if your file is not int16 use load/resample.m'])
skip = 0;
gains = [];
output = [filename '-gains'];

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+3) ' is not a property (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'output',
			output = varargin{i+1};
			if ~isstring_FMAT(output),
				error('Incorrect value for property ''output'' (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).');
			end
		case 'nchannels',
			nChannels = varargin{i+1};
			if ~isiscalar(nChannels,'>0'),
				error('Incorrect value for property ''nChannels'' (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).');
			end
		case 'gains',
			gains = varargin{i+1};
			if ~isdvector(gains,'>=0'),
				error('Incorrect value for property ''gains'' (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).');
			end
		case 'precision',
			precision = varargin{i+1};
			if ~isstring_FMAT(precision),
				error('Incorrect value for property ''precision'' (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).');
			end
		case 'skip',
			skip = varargin{i+1};
			if ~isiscalar(skip,'>=0'),
				error('Incorrect value for property ''skip'' (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).']);
	end
end

% Check output file
f = fopen(output,'r');
if f ~= -1,
	error(['Output file ''' output ''' exists.']);
end

% Check gains
gains = gains(:)';
nGains = length(gains);
if nGains == 1,
	gains = repmat(gains,1,nChannels);
	nGains = nChannels;
elseif nGains ~= nChannels,
	error('Incorrect number of gains (type ''help <a href="matlab:help ChangeBinaryGain">ChangeBinaryGain</a>'' for details).');
end

% Read in chunks, apply gains and save
offset = 0;
nSamples = floor(100000/nChannels);
while true,
	data = bz_LoadBinary(filename,'nChannels',nChannels,'offset',offset,'samples',nSamples,'precision',precision,'skip',skip);
	data = data .* repmat(gains,size(data,1),1);
	SaveBinary(output,data,'precision',precision,'mode','append');
	offset = offset + nSamples;
	if size(data,1) < nSamples, break; end
end
