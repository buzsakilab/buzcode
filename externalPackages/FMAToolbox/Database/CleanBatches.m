function CleanBatches

%CleanBatches - Delete completed batch jobs from memory.
%
% Batches cannot be removed from memory using 'clear' (they are based on Matlab
% timers, which are kept in memory even after the variables are cleared). This
% function will take the appropriate steps.
%
% This will not affect pending or running batch jobs.
%
%  USAGE
%
%    CleanBatches
%
%  SEE
%
%    See also StartBatch, GetBatch, BatchInfo, CancelBatch, Store, Recall.
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters

k = input('Make sure you have saved the results before you proceed!\nDo you wish to clean completed batch jobs (y/N)? ','s');
if isempty(k) || ~strcmp(lower(k(1)),'y'), return; end

delete(timerfindall('Tag','BatchJob','Running','off'));