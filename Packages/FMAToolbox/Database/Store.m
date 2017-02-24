function Store(variable,name,varargin)

%Store - Store variable in memory to avoid redundant computations.
%
%  This function is useful in the following scenario. When processing data
%  in batches, a function is called repeatedly with different parameters.
%  Because not all parameters change every time, some computations may be
%  repeated across successive function calls. For time-consuming computations,
%  (e.g. band-pass filtering of local field potentials) this represents a
%  significant waste of time.
%
%  Using Store and <a href="matlab:help Recall">Recall</a> solves this issue: you can 'store' the results of
%  a computation within a function, and 'recall' it next time the function
%  is called. Note that the variable is stored in memory, not on disk.
%
%  USAGE
%
%    Store(variable,name,key1,key2...)
%
%    variable       result to store
%    name           user-defined (arbitrary) name for the variable
%    key1,key2...   values of the parameters used to compute the variable
%
%  EXAMPLE
%
%    function y = ThetaPhaseLocking(tetrode,cluster,channel)
%
%      % First, try to recall the instantaneous phase for this
%      % channel (it may have been computed in a previous call)
%
%      phase = Recall('phase',channel);
%
%      if isempty(phase),
%
%        % Recalling failed, compute the phase now
%
%        lfp = GetLFP(channel);
%        theta = FilterLFP(lfp,'passband','theta');
%        phase = Phase(theta);
%
%        % Store the results for next time
%
%        Store(phase,'phase',channel);
%
%      end
%      ...
%
%  NOTE
%
%    Use parsimoniously to avoid saturating memory.
%
%  SEE
%
%    See also StartBatch, GetBatch, BatchInfo, CancelBatch, CleanBatches, Recall.
%

% Copyright (C) 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% (using weird names makes it unlikely that this global could be used somewhere else)
global storage_Store_Recall_981;

stack = dbstack;
if length(stack) < 2,
	functionName = 'CommandWindow';
else
	functionName = stack(2).name;
end

if ~exist('storage_Store_Recall_981'),
	storage_Store_Recall_981 = struct;
end

% Remove previously stored variable with this name, if any
if isfield(storage_Store_Recall_981,functionName),
	topFieldName = ['storage_Store_Recall_981.' functionName];
	topField = eval(topFieldName);
	if isfield(topField,name), eval([topFieldName ' = rmfield(' topFieldName ',''' name ''');']); end
end

% Store variable
storage_Store_Recall_981 = setfield(storage_Store_Recall_981,functionName,name,'x',variable);

% Store keys
for i = 1:length(varargin),
	storage_Store_Recall_981 = setfield(storage_Store_Recall_981,functionName,name,['key' int2str(i)],varargin{i});
end
