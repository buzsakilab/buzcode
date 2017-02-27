function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: vertcat.m 8776 2013-11-14 09:04:48Z roboos $


o = cat(1,varargin{:});
return;
