function tf = isfield(this,field)
% Isfield method for GIfTI objects
% FORMAT tf = isfield(this,field)
% this   -  GIfTI object
% field  -  string of cell array
% tf     -  logical array
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: isfield.m 8776 2013-11-14 09:04:48Z roboos $

tf = ismember(field, fieldnames(this));