function CancelBatch(batch)

%CancelBatch - Cancel batch job.
%
%  USAGE
%
%    CancelBatch(batch)
%
%    batch          batch parameter returned by <a href="matlab:help StartBatch">StartBatch</a>
%
%  SEE
%
%    See also StartBatch, GetBatch, BatchInfo, CleanBatches, Store, Recall.
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help CancelBatch">CancelBatch</a>'' for details).');
end
if prod(size(batch)) ~=1 || ~isa(batch,'timer') || ~isvalid(batch) || ~strcmp(get(batch,'Tag'),'BatchJob'),
	error('Incorrect parameter (type ''help <a href="matlab:help CancelBatch">CancelBatch</a>'' for details).');
end

k = input('Are you sure you wish to cancel this batch (y/N)? ','s');
if isempty(k) || ~strcmp(lower(k(1)),'y'), return; end

delete(batch);