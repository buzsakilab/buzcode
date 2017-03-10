function batch = StartBatch(mfile,bfile,varargin)

%StartBatch - Start a new batch job.
%
% Batch processing is useful if you need to repeatedly perform a given analysis
% on different data sets. Using batch processing is easier and more reliable
% than writing a simple 'for' loop:
%
%  * You must only supply the code to run the analysis on one data set,
%    and a text file listing the data sets, including specific parameters
%
%  * Running the batch on a subset of the data is as easy as commenting
%    out the appropriate lines in this file
%
%  * Batch processing will not be interrupted when the analysis of one
%    particular data set generates an error
%
%  * Errors can be logged to disk
%
%  * Processing can be delayed to a later time (e.g. during the night)
%
% Results can either be kept in memory using the batch output variable, or
% stored in a central database. The latter is the recommended solution, as it
% is far more powerful and reliable. See EXAMPLES below.
%
% During batch processing, figures are typically hidden so they do not keep
% popping up as the user is performing other tasks as the batch is running in
% the background (although this is configurable). In order to keep figures
% hidden, batch functions should not use figure(h) or axes(a) which make the
% figure visible. Two helper functions, scf(h) and sca(h), can be used instead
% to set the current figure or axes without making the figure visible.
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
%     'log'         error log file (default = none) (see NOTE below)
%    =========================================================================
%
%  NOTE
%
%    Log file names can include wildcards to indicate current date and time:
%
%      %y    year
%      %m    month
%      %d    day
%      %t    time
%
%  EXAMPLE 1 - RESULTS KEPT IN BATCH OUTPUT VARIABLE
%
%    This is the simpler solution, but it is not as reliable (e.g. the results
%    can be lost if Matlab crashes) and powerful (e.g. the results cannot be
%    shared between users) as using a database. See EXAMPLE 2 below.
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
%        lfp = GetLFP(channel);
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
%    To log errors:
%
%        b = StartBatch(@CountRipples,'batch.txt','log','Errors-%y%m%d.log');
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
%  EXAMPLE 2 - RESULTS STORED IN CENTRAL DATABASE
%
%    Let us rewrite the previous example using a central database. First, we need
%    to create the database:
%
%        % This assumes your login/password are stored
%        % in a configuration file in your home directory
%        DBConnect;
%        DBCreate('Ripples_%y%m%d');
%
%    The function used to process the data would now look like:
%
%        function CountRipples(session,channel)
%
%        SetCurrentSession(session);
%        lfp = GetLFP(channel);
%        fil = FilterLFP(lfp,'bandpass','ripples');
%        ripples = FindRipples(fil);
%        n = size(ripples,1);
%        eid = [session '-' channel];
%        DBAddVariable(n,eid,'Ripple Count','','',{'CountRipples'});
%
%    This would store the number of ripples for each session, as well as useful
%    information such as who stored these data, when, using which function (the
%    actual code, not just the function name, is stored in the database).
%
%    To get the results:
%
%        n = DBGetValues;
%
%  SEE
%
%    See also Database, GetBatch, BatchInfo, CancelBatch, CleanBatches, Store,
%    Recall, Hide, sca, scf.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
delay = 0;
hide = 'off';
log = '';

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
    case 'log',
      log = varargin{i+1};
      if ~isstring_FMAT(log),
        error('Incorrect value for property ''log'' (type ''help <a href="matlab:help StartBatch">StartBatch</a>'' for details).');
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

% Parse batch file
b = ParseBatch(bfile);

% Open log file
b.log = -1;
if ~isempty(log),
	% Insert date
	log = InsertDate(log);
	b.log = fopen(log,'w');
	if b.log == -1, error(['Could not open log file ''' log '''.']); end
end

% Set batch function
[~,mfile] = fileparts(mfileName);
b.mfile = mfile;

b.hide = strcmp(hide,'on');
if b.hide,
	warning('Figures will be hidden (type ''help <a href="matlab:help Hide">Hide</a>'' for details).');
end

% Start timer
batch = timer('TimerFcn',{@RunBatch,b},'StartDelay',delay*60,'Tag','BatchJob');
start(batch);
