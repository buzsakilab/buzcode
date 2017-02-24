function MClustCutterUndoRedo()

% MClustCutterUndoRedo()
%
% ADR 2008
%
% pop the redo stack
%
% Part of the new undo mechanism in MClust 3.5.
% Modified from code by JCJ 2007
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5

global MClust_Undo MClust_Redo
global MClust_Clusters  MClust_Colors 

% Increment Undo pointer if further undo information exists and last non-undid state exists.
if isempty(MClust_Redo)
    warndlg('No steps to redo.')
    return
end

% push the undo stack
U = MClustCutterUndoGetCurrentState(MClust_Redo{end}.funcname);
MClust_Undo{end+1} = U;

% pop redo stack
MC_Redo = MClust_Redo{end}; 
MClust_Redo(end) = [];
MClust_Clusters = MC_Redo.clusters;
MClust_Colors = MC_Redo.colors; 

% function name redone is store under current undo
msgbox(['Re-applied function ', MC_Redo.funcname], 'Redo', 'none', 'non-modal'); % made non-modal to remove error sounds

global MClust_FilesWrittenYN
MClust_FilesWrittenYN = 'no';

% update tooltips
MClustCutterUndoUpdateTooltip();

% count down autosave
MClustCutterStepAutosave;
