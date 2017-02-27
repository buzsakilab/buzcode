function hdr = empty_hdr
% Create an empty NIFTI-1 header
% FORMAT hdr = empty_hdr
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: empty_hdr.m 8776 2013-11-14 09:04:48Z roboos $


org = niftistruc;
for i=1:length(org),
    hdr.(org(i).label) = org(i).def;
end;

