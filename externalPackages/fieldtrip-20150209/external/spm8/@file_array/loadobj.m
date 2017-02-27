function b = loadobj(a)
% loadobj for file_array class
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: loadobj.m 8776 2013-11-14 09:04:48Z roboos $

if isa(a,'file_array')
    b = a;
else
    a = permission(a, 'rw');
    b = file_array(a);
end