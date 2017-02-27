function out = double(fa)
% Convert to double precision
% FORMAT double(fa)
% fa - a file_array
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: double.m 8776 2013-11-14 09:04:48Z roboos $

out = double(numeric(fa));

