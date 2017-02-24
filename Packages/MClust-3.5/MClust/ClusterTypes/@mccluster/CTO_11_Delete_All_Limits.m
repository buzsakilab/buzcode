function [MCC, redraw, rekey, undoable] = CTO_Delete_AllLimits(MCC, varargin)

% [MCC, redraw, rekey, undoable] = CTO_Delete_AllLimits(MCC, varargin)
% REQUIRES
%
% PL 2000-2002
% ADR 2008
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

global MClust_FeatureSources

redraw= true; rekey = false; undoable = true;

MCC.xdimNames = {}; 
MCC.ydimNames = {}; 
MCC.xdimSources = {};
MCC.ydimSources = {};
MCC.cx = {}; 
MCC.cy = {};

% recalculate myPoints
MCC.recalc = 1; % things have changed;  
MCC.myPoints = FindInCluster(MCC);


