function info = BatchInfo(batch)

%BatchInfo - Get batch job information.
%
%  USAGE
%
%    info = BatchInfo(batch)
%
%    batch          batch parameter returned by <a href="matlab:help StartBatch">StartBatch</a>
%
%  OUTPUT
%
%    info.mfile     batch function
%    info.data      batch data
%    info.done      whether the job has completed
%
%  SEE
%
%    See also StartBatch, GetBatch, CancelBatch, CleanBatches, Store, Recall.
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help BatchInfo">BatchInfo</a>'' for details).');
end
if prod(size(batch)) ~=1 || ~isa(batch,'timer') || ~isvalid(batch) || ~strcmp(get(batch,'Tag'),'BatchJob'),
	error('Invalid batch job (type ''help <a href="matlab:help BatchInfo">BatchInfo</a>'' for details).');
end

job = get(batch,'TimerFcn');
info.mfile = job{2}.mfile;
info.data = job{2}.field;
info.done = strcmp(lower(get(batch,'Running')),'off');

