function [MCC, redraw, rekey, undoable] = CTO_02_Delete_Limit(MCC, varargin)

% MCCluster/DeleteLimit(MCC, varargin)
%
% REQUIRES
%    xdim - xdimension
%    ydim - ydimension
%
% PL 2000-2002
% ADR 2008
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

global MClust_CurrentFeatureNames

extract_varargin;

redraw= true; rekey = false; undoable = true;

if isempty(MCC.xdimNames) || isempty(MCC.ydimNames) % ADR 2008
    return
end

fx = strmatch(MClust_CurrentFeatureNames{1}, MCC.xdimNames);
fy = strmatch(MClust_CurrentFeatureNames{2}, MCC.ydimNames);
iL = intersect(fx,fy);

if ~isempty(iL)
    MCC.xdimNames(iL) = [];
    MCC.ydimNames(iL) = [];
    MCC.xdimSources(iL,:) = [];
    MCC.ydimSources(iL,:) = [];
    MCC.cx(iL) = [];
    MCC.cy(iL) = [];
end

