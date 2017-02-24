function [MCC, redraw, rekey, undoable] = DeletePoints(MCC, varargin)

% [MCC, redraw, rekey, undoable] = CTO_Delete_AllPoints(MCC, varargin)
%
% INPUTS
%     MCC - a MCCluster
% needs xdim, ydim, drawingAxis from parent
%
% OUTPUTS
%     MCC - The updated cluster
%
% 
% ncst 26 Nov 02
% ADR 2008
%

redraw= true; rekey = false; undoable = true;

global MClust_CurrentFeatureData

extract_varargin;

if isempty(drawingAxis)
	errordlg('No drawing axis available.', 'Error');
	return
end

axes(drawingAxis);
[chx,chy] = DrawConvexHull;

if ~isempty(chx) && ~isempty(chy)
    f = MCC.myOrigPoints;

    f_del = find(InPolygon(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), chx, chy));
    MCC.myOrigPoints(f_del) = [];

    MCC.recalc = 1;
end