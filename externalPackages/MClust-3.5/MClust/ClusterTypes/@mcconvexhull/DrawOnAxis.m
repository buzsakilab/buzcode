function DrawConvexHull(MCC, xdim, ydim, color, axisHandle)

% mcconvexhull/DrawConvexHull(MCC, xd, yd, color, axisHandle)
%
% ADR 1998-2008
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.5.

global MClust_CurrentFeatureNames

fx = strmatch(MClust_CurrentFeatureNames{1}, MCC.xdimNames);
fy = strmatch(MClust_CurrentFeatureNames{2}, MCC.ydimNames);
iL = intersect(fx,fy);
if ~isempty(iL)
    axes(axisHandle)
    h = plot(MCC.cx{iL}, MCC.cy{iL}, '-');
    set(h, 'Color', color);
end
