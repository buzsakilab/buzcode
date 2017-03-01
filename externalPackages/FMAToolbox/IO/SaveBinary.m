function SaveBinary(filename,data,varargin)

%SaveBinary - Save data in a multiplexed binary file.
%
%  USAGE
%
%    SaveBinary(filename,data,<options>)
%
%    filename       file to save to
%    data           data to save
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'precision'   sample precision (default = 'int16')
%     'mode'        either 'new' to create a new file (default) or 'append'
%                   to add to an existing file
%    =========================================================================

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
% precision = 'int16'; % THIS SHOULD NEVER BE ASSUMED
precision = class(data);
mode = 'new';

if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SaveBinary">SaveBinary</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help <a href="matlab:help SaveBinary">SaveBinary</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'precision',
      precision = varargin{i+1};
      if ~isstring_FMAT(precision),
        error('Incorrect value for property ''precision'' (type ''help <a href="matlab:help SaveBinary">SaveBinary</a>'' for details).');
      end
    case 'mode',
      mode = lower(varargin{i+1});
      if ~isstring_FMAT(mode,'new','append'),
        error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help SaveBinary">SaveBinary</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SaveBinary">SaveBinary</a>'' for details).']);
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

if strcmp(mode,'new'),
	f = fopen(filename,'w');
	status = fseek(f,0,'bof');
else
	f = fopen(filename,'a');
	status = fseek(f,0,'eof');
end
if status ~= 0,
	fclose(f);
	error('Could not start writing (possible reasons include insufficient access permissions).');
end

fwrite(f,data',precision);

fclose(f);
