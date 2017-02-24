function currentStateForUndo = MClustCutterUndoGetCurrentState(funcname)

% currentStateForUndo = MClustCutterUndoGetCurrentState(funcname)
%
% ADR 2008
%
% Part of the new undo mechanism in MClust 3.5.
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

global MClust_Colors
global MClust_Clusters

currentStateForUndo.colors = MClust_Colors; 
currentStateForUndo.clusters = MClust_Clusters;
currentStateForUndo.funcname = funcname;
