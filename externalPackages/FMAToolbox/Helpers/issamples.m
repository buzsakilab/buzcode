%issamples - Test if parameter is a list of samples satisfying an optional list of tests.
%
%  A list of samples is a vector or matrix of doubles containing at least one column
%  for timestamps, and optionally more for continuous variables.
%
%  USAGE
%
%    test = issamples(x,test1,test2,...)
%
%    x              parameter to test
%    test1...       optional list of additional tests (see examples below)
%
%  EXAMPLES
%
%    % Test if x is a list of samples
%    issamples(x)
%
%    % Special test: test if x is a point process (i.e. 1 column)
%    issamples(x,'#0')
%
%    % Special test: test if x contains 3 variables (i.e. 4 columns)
%    issamples(x,'#3')
%
%  SEE ALSO
%
%    See also isdmatrix, isdvector, isdscalar, isimatrix, isivector, isiscalar,
%    isstring_FMAT, islscalar, islvector, islmatrix.
%

% Copyright (C) 2010-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = issamples(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help issamples">issamples</a>'' for details).');
end

% Test: double, vector
test = isa(x,'double');

% Optional tests
for i = 1:length(varargin),
	try
		if varargin{i}(1) == '#',
			if size(x,2)-1 ~= str2num(varargin{i}(2:end)), test = false; return; end
		end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help issamples">issamples</a>'' for details).']);
	end
end
