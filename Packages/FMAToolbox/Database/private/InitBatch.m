function b = InitBatch(bfile)

%InitBatch - Initialize a new batch job.
%
%  USAGE
%
%    b = InitBatch(bfile)
%
%    bfile          batch file listing the parameters for each iteration
%
%  SEE
%
%    See also StartBatch, ShowBatch
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help InitBatch">InitBatch</a>'' for details).');
end

% Make sure batch file is valid
s = dbstack;
if ~isstring_FMAT(bfile) || ~exist(bfile,'file'),
	error(['Batch file not found (type ''help <a href="matlab:help ' s.name '">' s.name '</a>'' for details).']);
end

% Open batch file
f = fopen(bfile,'r');
if f == -1, error(['Could not open file ''' bfile '''.']); end

item = 1;
l = 1;
nFields = 0;
% Read batch file
while ~feof(f),
	field = 0;
	% Read next line
	line = fgetl(f);
	while ~isempty(line),
		% Get next field
		[token,line] = strtok(line);
		% Skip rest of line if this is a comment
		if isempty(token) | token(1) == '%', break; end
		field = field + 1;
		% Determine if this is a number, vector or string
		n = str2num(token);
		if ~isempty(n),
			% It is a number, convert it to numerical format
			b.field{item,field} = n;
		elseif token(1) ~= '[',
			% It is a string, keep it as it is
			b.field{item,field} = token;
		else
			% It is a vector (or matrix)
			[token,line] = strtok([token line],'[]');
			n = str2num(token);
			if isempty(n),
				error(['Incorrect vector or matrix in batch file (line ' int2str(l) ', element ' int2str(field) ').']);
			end
			b.field{item,field} = n;
			line = line(2:end);
		end
	end
	if field > 0,
		% Make sure all lines have the same number of elements
		if item == 1,
			nFields = field;
		elseif field ~= nFields,
			if l == 2,
				error(['Incoherent numbers of elements in batch file (' int2str(nFields)  ' elements on line 1, but ' int2str(field) ' elements on line 2).']);
			else
				error(['Incoherent numbers of elements in batch file (' int2str(nFields)  ' elements on lines 1-' int2str(l-1) ', but ' int2str(field) ' elements on line ' int2str(l) ').']);
			end
		end
		% Next item
		b.line(item) = l;
		item = item + 1;
	end
	% Update line counter (for error messages)
	l = l + 1;
end

% Close file
fclose(f);

% Reset iterator
b.currentItem = 0;
b.currentField = 0;

% Empty fields
b.log = -1;
b.mfile = [];
