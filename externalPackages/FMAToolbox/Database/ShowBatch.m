function b = ShowBatch(bfile)

%ShowBatch - Show data sets in a batch job.
%
%  USAGE
%
%    b = ShowBatch(bfile)
%
%    bfile          batch file listing the parameters for each iteration
%
%    See also StartBatch, GetBatch, BatchInfo, CancelBatch, CleanBatches, Store, Recall.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ShowBatch">ShowBatch</a>'' for details).');
end

% Make sure the batch file exists
if ~isstring(bfile) || ~exist(bfile,'file'),
	error('Batch file not found (type ''help <a href="matlab:help ShowBatch">ShowBatch</a>'' for details).');
end

% Open batch file
f = fopen(bfile,'r');
if f == -1, error(['Could not open file ''' bfile '''.']); end

b = ParseBatch(bfile);
b.field

