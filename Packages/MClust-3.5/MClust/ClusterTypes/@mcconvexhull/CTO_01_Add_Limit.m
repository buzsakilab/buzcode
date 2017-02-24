function [MCC, redraw, rekey, undoable] = AddLimit(MCC, varargin)

% MCCluster/AddLimit(MCC, varargin)
%
% REQUIRES
%    xdim - xdimension
%    ydim - ydimension
%    axisHandle - axis on which to draw convex hull
%
% ADR 1998
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

global MClust_CurrentFeatureNames MClust_FeatureSources
extract_varargin;

redraw= true; rekey = false; undoable = true;

if isempty(drawingAxis)
	errordlg('No drawing axis available.', 'Error');
	return
end

fx = strmatch(MClust_CurrentFeatureNames{1}, MCC.xdimNames);
fy = strmatch(MClust_CurrentFeatureNames{2}, MCC.ydimNames);
iL = intersect(fx,fy);

if isempty(iL)
   iL = length(MCC.xdimNames)+1;
end

axes(drawingAxis);
[chx,chy] = DrawConvexHull;

if ~isempty(chx) && ~isempty(chy)
    MCC.xdimNames{iL} = MClust_CurrentFeatureNames{1};
    MCC.xdimSources{iL,1} = MClust_FeatureSources{xdim,1};
    MCC.xdimSources{iL,2} = MClust_FeatureSources{xdim,2};
    MCC.ydimNames{iL} = MClust_CurrentFeatureNames{2};
    MCC.ydimSources{iL,1} = MClust_FeatureSources{ydim,1};
    MCC.ydimSources{iL,2} = MClust_FeatureSources{ydim,2};
    MCC.cx{iL} = chx;
    MCC.cy{iL} = chy;
end

