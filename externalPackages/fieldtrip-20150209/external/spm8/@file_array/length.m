function l = length(x)
% Overloaded length function for file_array objects
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: length.m 8776 2013-11-14 09:04:48Z roboos $


l = max(size(x));

