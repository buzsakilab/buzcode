function [redraw, rekey, undoable] = ShowContourPlot

% [redraw, rekey, undoable] = ShowContourPlot
%
% INPUTS
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
%
% Changes type of iClust
%
% ADR 2008
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

redraw = false; rekey = false; undoable = false;

contourWindow = findobj('Type', 'figure', 'Tag', 'ContourWindow');
drawingFigHandle = findobj('Type', 'figure', 'Tag', 'CHDrawingAxisWindow');
if isempty(contourWindow)
	contourWindow = figure('NumberTitle', 'off', 'Name', 'Contour Plot', 'Tag', 'ContourWindow');
end
mkContours(drawingFigHandle, 'figHandle', contourWindow);