function [MCC, redraw, rekey, undoable] = AddPoints(MCC, varargin)

% f = AddPoints(MCC, varargin)
%
% INPUTS
%     MCC - a MCCluster
% needs xdim, ydim, drawingAxis from parent
%
% c
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
    f = MCC.myPoints;
    f_new = find(InPolygon(MClust_CurrentFeatureData(:,1), MClust_CurrentFeatureData(:,2), chx, chy));	 
    MCC.myPoints = unique([f; f_new]);
end