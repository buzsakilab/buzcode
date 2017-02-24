%RunBatch - Run a batch job. This should *not* be called directly.
%
% Run a batch job. This function is called automatically by the
% batch timer upon expiration of the required delay.
%
%  USAGE
%
%    RunBatch(timer,event,b)
%
%    timer          Matlab timer object
%    event          Matlab timer event type
%    b              batch structure
%

% Copyright (C) 2007-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function RunBatch(timer,event,b)

if b.hide,
	Hide('on');
	warning('off','FMAToolbox:Hide:FigureHidden');
end

while true,
	% Get next item
	[b,item] = GetNextItem(b);
	if isempty(item), break; end
	% Start building command line
	i = 1;
	clear('args');
	while true,
		% Get next field
		[b,field] = GetNextField(b);
		if isempty(field), break; end
		args{i} = field;
		i = i + 1;
	end

	try
		% Call batch function
		outputs = cell(1,nargout(b.mfile));
		[outputs{:}] = feval(b.mfile,args{:});
		% Append results to 'UserData' field of timer object (this is a cell array)
		data = timer.UserData;
		if isempty(data),
			data = {outputs{:}};
		else
			data(end+1,:) = {outputs{:}};
		end
		timer.UserData = data;
	catch
		fprintf(2,['Batch: error processing item ' int2str(item) '\n']);
		e = lasterror;
		fprintf(2,[' ' e.message '\n']);
		fprintf(2,[' Error in ==> ' e.stack(1).name ' at ' int2str(e.stack(1).line) '\n']);
	end
end

if b.hide,
	Hide('off');
	Hide('none');
	warning('on','FMAToolbox:Hide:FigureHidden');
end

