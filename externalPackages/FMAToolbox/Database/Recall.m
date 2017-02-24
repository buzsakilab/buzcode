function variable = Recall(name,varargin)

%Recall - Recall variable from memory to avoid redundant computations.
%
%  This function is useful in the following scenario. When processing data
%  in batches, a function is called repeatedly with different parameters.
%  Because not all parameters change every time, some computations may be
%  repeated across successive function calls. For time-consuming computations,
%  (e.g. band-pass filtering of local field potentials) this represents a
%  significant waste of time.
%
%  Using <a href="matlab:help Store">Store</a> and Recall solves this issue.
%
%  USAGE
%
%    x = Recall(name,key1,key2...)
%
%    name           user-defined (arbitrary) name for the variable
%    key1,key2...   values of the parameters used to compute the variable
%
%  EXAMPLE
%
%    function y = ThetaPhaseLocking(tetrode,cluster,channel)
%      phase = Recall('phase',channel);
%      if isempty(phase),
%        lfp = GetLFP(channel);
%        theta = FilterLFP(lfp,'passband','theta');
%        phase = Phase(theta);
%        Store(phase,'phase',channel);
%      end
%      ...
%
%  NOTE
%
%    Use parsimoniously to avoid saturating memory.
%
%  SEE
%
%    See also StartBatch, GetBatch, BatchInfo, CancelBatch, CleanBatches, Store.
%

% Copyright (C) 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% (using weird names makes it unlikely that this global could be used somewhere else)
global storage_Store_Recall_981;

variable = [];

stack = dbstack;
if length(stack) < 2,
	functionName = 'CommandWindow';
else
	functionName = stack(2).name;
end

try
	for i = 1:length(varargin),
		x = getfield(storage_Store_Recall_981,functionName,name,['key' int2str(i)]);
		if x ~= varargin{i}, return; end
	end
	variable = getfield(storage_Store_Recall_981,functionName,name,'x');
catch
	return
end
