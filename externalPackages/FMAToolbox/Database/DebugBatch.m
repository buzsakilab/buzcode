function DebugBatch(mfile,bfile,item,varargin)

%DebugBatch - Help debug a batch job.
%
%  This function can be used to help debug a batch. It retrieves the parameters
%  corresponding to the required item from the batch file, and either starts
%  the batch function with these parameters in debug mode, or simply assigns
%  these parameters to new variables in matlab's workspace so they can be
%  directly manipulated.
%
%  Set breakpoints in the batch function before running DebugBatch. If no
%  breakpoints have been defined, DebugBatch will stop at the first line.
%
%  USAGE
%
%    DebugBatch(mfile,bfile,item,<options>)
%
%    mfile          batch function (M-file name or function handle)
%    bfile          batch file listing the parameters for each iteration
%    item           item number in batch function
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'        run the batch function ('run', default) or simply set
%                   variables in matlab's workspace ('set')
%    =========================================================================
%
%  OUTPUT
%
%    Instantiates the variables of the batch function using the values in the
%    batch file ('set' mode).
%
%  SEE
%
%    See also StartBatch.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
mode = 'run';

% Check number of parameters
if nargin < 3,
	error(['Incorrect number of parameters (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
end

% Batch function name
if isa(mfile,'function_handle'),
	mfile = func2str(mfile);
end

% Check batch file and function are valid
if ~isstring(bfile) || ~exist(bfile,'file'),
	error(['Batch file not found (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
end
if isempty(which(mfile)),
	error(['Batch function not found (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
end

% Item
if ~isiscalar(item),
	error(['Incorrect item (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'mode',
      mode = lower(varargin{i+1});
      if ~isstring(mode,'run','set'),
        error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
  end
end

% Open batch file
f = fopen(bfile,'r');
if f == -1, error(['Could not open file ''' bfile '''.']); end

% Parse batch file
b = ParseBatch(bfile);

if strcmp(mode,'run'),

	% Check breakpoints
	if length(dbstatus(mfile)) == 0,
		eval(['dbstop in ' mfile ' at 1;']);
	end
	
	% Run batch function with appropriate parameters
	for i = 1:size(b.field,2),
		args{i} = b.field{item,i};
	end
	outputs = cell(1,nargout(mfile));
	[outputs{:}] = feval(mfile,args{:})

else

	% Open batch function
	if isa(mfile,'function_handle'),
		mfile = func2str(mfile);
	end
	f = fopen(which(mfile));

	% Find first line containing the (uncommented) 'function' keyword
	found = [];
	while isempty(found),
		line = fgets(f);
		if line == -1, error('Could not find function definition in batch function.'); end
		line = regexprep(line,'%.*function.*','');
		found = regexp(line,'.*function[^(]*\(([^)]*)\).*');
	end
	fclose(f);

	% Extract parameter names, and assign them in 'base' workspace
	parameters = regexprep(line,'.*function[^(]*\(([^)]*)\).*','$1');
	parameters = regexp(parameters,'[^,]*','match');
	for i = 1:length(parameters),
		assignin('base',parameters{i},b.field{item,i});
	end
	
end