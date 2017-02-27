function out = ndims(fa)
% Number of dimensions
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: ndims.m 8776 2013-11-14 09:04:48Z roboos $


out = size(fa);
out = length(out);

