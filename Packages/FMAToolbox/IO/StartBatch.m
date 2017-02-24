function batch = StartBatch(mfile,bfile,varargin)

%StartBatch - Start a new batch job.
%
% Batch jobs are useful if you need to repeatedly run a given function on
% different parameters, but a simple 'for' loop would not be adequate.
% Using batch jobs, the parameters for each iteration are read from a text
% file, called a 'batch file'. See Example below.
%
%  USAGE
%
%    batch = StartBatch(mfile,bfile,<options>)
%
%    mfile          batch function (M-file name or function handle)
%    bfile          batch file listing the parameters for each iteration
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'delay'       delay (in minutes) before execution starts
%     'hide'        by default, new figures are hidden to speed up batch
%                   processing (use 'off' to cancel this behavior)
%    =========================================================================
%
%  EXAMPLE
%
%    Suppose you wish to count the number of 'ripples' in a set of sleep
%    sessions recorded from different rats. You could create the file
%    'batch.txt', where you would list the session files to process, as
%    well as the LFP channel to use in each case to detect ripples:
%
%        Rat-01-20100915-01-sleep     5
%        Rat-01-20100920-01-sleep     5
%        Rat-02-20101001-01-sleep     3
%        Rat-01-20101005-04-rest      1
%        ...
%
%    The function used to process the data would look like:
%
%        function n = CountRipples(session,channel)
%
%        SetCurrentSession(session);
%        lfp = StartBatch(channel);
%        fil = FilterLFP(lfp,'bandpass','ripples');
%        ripples = FindRipples(fil);
%        n = size(ripples,1);
%
%    To start processing the data in one hour, you would then call:
%
%        b = StartBatch(@CountRipples,'batch.txt','delay',60);
%
%    To get the results:
%
%        n = GetBatch(b);
%
%  NOTE
%
%    For technical reasons, the batch function must have a fixed number of output
%    parameters, e.g. it cannot have 'varargout' as one of its output parameters.
%    It is possible to circumvent this limitation if you know in advance the number
%    of output parameters that the batch function would actually return in your
%    particular case. Assuming this function is called 'VariableProcess', takes two
%    inputs and will always return three outputs for your data, you could define:
%
%        function [x,y,z] = FixedProcess(u,v)
%
%        [x,y,z] = VariableProcess(u,v);
%
%    and then use this as the batch function.
%
%  SEE
%
%    See also GetBatch, BatchInfo, CancelBatch, CleanBatches, Store, Recall.
%

% Copyright (C) 2007-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
delay = 0;
hide = 'on';

% Check number of parameters
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'delay',
      delay = varargin{i+1};
      if ~isiscalar(delay,'>0'),
        error('Incorrect value for property ''delay'' (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details).');
      end
    case 'hide',
      hide = lower(varargin{i+1});
      if ~isstring_FMAT(hide,'on','off'),
        error('Incorrect value for property ''hide'' (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details).']);
  end
end

% Compatibility with previous versions: reverse parameter order if necessary
if isa(bfile,'function_handle') || (isstring_FMAT(mfile) && isempty(which(mfile))),
	tmp = bfile;
	bfile = mfile;
	mfile = tmp;
end

% Batch function name
if isa(mfile,'function_handle'),
	mfileName = func2str(mfile);
else
	mfileName = mfile;
end

% Check batch file and function are valid
if ~isstring_FMAT(bfile) || ~exist(bfile,'file'),
	error('Batch file not found (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details).');
end
if isempty(which(mfileName)),
	error('Batch function not found (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details).');
end

% Make sure M-file has a defined number of output parameters
if nargout(mfile) == -1,
	error(['''' mfileName ''' has a variable number of output arguments (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details). ']);
end

% Open batch file
f = fopen(bfile,'r');
if f == -1, error(['Could not open file ''' bfile '''.']); end

item = 1;
% Read batch file
while ~feof(f),
	field = 0;
	% Read next line
	line = fgetl(f);
	while ~isempty(line),
		% Get next field
		[token,line] = strtok(line);
		% Skip rest of line if this is a comment
		if isempty(token) | token(1) == '%', break; end
		field = field + 1;
		% Determine if this is a number
		n = str2num(token);
		if isempty(n),
			% It is a string, keep it as it is
			b.field{item,field} = token;
		else
			% It is a number, convert it to numerical format
			b.field{item,field} = n;
		end
	end
	if field > 0, item = item + 1; end
end

% Close file
fclose(f);

% Reset iterator
b.currentItem = 0;
b.currentField = 0;

% Set batch function
[unused,mfile] = fileparts(mfileName);
b.mfile = mfile;

b.hide = strcmp(hide,'on');
if b.hide,
	warning('Figures will be hidden (type ''help <a href="matlab:help Hide">Hide</a>'' for details).');
end

% Start timer
batch = timer('TimerFcn',{@RunBatch,b},'StartDelay',delay*60,'Tag','BatchJob');
start(batch);
