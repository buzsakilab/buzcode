function chkinputdatatype(varargin)
%CHKINPUTDATATYPE Check that all inputs are double

%   Copyright 2009-2013 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

for n = 1:nargin
    if ~isa(varargin{n},'double')
        error(message('signal:chkinputdatatype:NotSupported'));
    end
end



% [EOF]
