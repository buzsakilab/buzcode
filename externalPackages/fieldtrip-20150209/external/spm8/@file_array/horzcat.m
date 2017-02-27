function o = horzcat(varargin)
% Horizontal concatenation of file_array objects
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: horzcat.m 8776 2013-11-14 09:04:48Z roboos $

o    = cat(2,varargin{:});
return;

