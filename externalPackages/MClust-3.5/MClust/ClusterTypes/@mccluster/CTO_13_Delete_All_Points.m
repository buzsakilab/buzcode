function [MCC, redraw, rekey, undoable] = CTO_Delete_AllPoints(MCC, varargin)

% [MCC, redraw, rekey, undoable] = CTO_Delete_AllPoints(MCC, varargin)
%
% REQUIRES
%
% ADR 2008
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

extract_varargin;

redraw= true; rekey = false; undoable = true;

MCC.myPoints = [];
MCC.myOrigPoints = [];

% recalculate myPoints
MCC.recalc = 1; % things have changed;  


