function [x,varargout] = ZeroToOne(x,varargin)

%ZeroToOne - Normalize values in [0,1].
%
%  USAGE
%
%    [y,b0,b1...] = ZeroToOne(x,a0,a1,...)
%
%    x              array to normalize (vector OR matrix)
%    a0...          additional inputs to transform using the same
%                   scale as x

% Copyright (C) 2008-2011 by Michaël Zugaro, modified GG 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ZeroToOne">ZeroToOne</a>'' for details).');
end

if nargin ~= nargout,
	error('Different numbers of input and output parameters (type ''help <a href="matlab:help ZeroToOne">ZeroToOne</a>'' for details).');
end

if isvector(x)
    m = min(x);
    M = max(x);
    x = (x-m)/(M-m);
elseif ismatrix(x)
    m = min(min(x));
    M = max(max(x));
    x = (x-m)./(M-m);
end

for i = 1:nargin-1,
	varargout{i} = (varargin{i}-m)/(M-m);
end
