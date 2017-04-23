function [type] = bz_IsEventDatatype(varargin)

% this function checks that the input variable is of the buzcode event
% datatype format

% David Tingley, 2017

data = varargin{1};

if isstruct(data) & isfield(data,'times')
    type = 1;
else
    type = 0;
end



    
    