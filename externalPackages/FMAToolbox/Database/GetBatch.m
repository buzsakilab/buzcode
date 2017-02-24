function varargout = GetBatch(batch)

%GetBatch - Get batch job output.
%
%  USAGE
%
%    [output1,output2,...] = GetBatch(batch)
%
%    batch          batch parameter returned by <a href="matlab:help StartBatch">StartBatch</a>
%
%  OUTPUT
%
%    output1...     one matrix or cell array for each output parameter of the
%                   batch function, listing the result of each individual iteration
%                   on a separate line
%
%  EXAMPLE
%
%    % Start a batch in 4 hours
%    >> b = StartBatch(@SomeComputation,'batchfile.txt',4*60);
%
%    % Come back in 3 hours and try getting the results
%    >> output = GetBatch(b);
%    Warning: Batch not finished...
%
%    % Come back the following morning and get the results
%    >> output = GetBatch(b);
%
%  SEE
%
%    See also StartBatch, BatchInfo, CancelBatch, CleanBatches, Store, Recall.
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help GetBatch">GetBatch</a>'' for details).');
end
if prod(size(batch)) ~=1 || ~isa(batch,'timer') || ~isvalid(batch) || ~strcmp(get(batch,'Tag'),'BatchJob'),
	error('Invalid batch job (type ''help <a href="matlab:help GetBatch">GetBatch</a>'' for details).');
end

if strcmp(get(batch,'Running'),'on'),
	warning('Batch not finished...');
	output = [];
else
	output = get(batch,'UserData');
end

if nargout > size(output,2),
	error('Too many output parameters (type ''help <a href="matlab:help GetBatch">GetBatch</a>'' for details).');
end

for i = 1:max([1 nargout]),
	if isnumeric(output{1,i}) || islogical(output{1,i}),
		% If this parameter is numeric or logical and all instances have the same size, concatenate,
		% otherwise return a cell array
		if any(diff(cellfun('size',output,2))),
			varargout{i} = {output{:,i}}';
		else
			varargout{i} = vertcat(output{:,i});
		end
	else
		varargout{i} = {output{:,i}}';
	end
end

