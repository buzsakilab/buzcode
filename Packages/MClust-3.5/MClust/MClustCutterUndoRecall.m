function MClustCutterUndoRecall()

% MClustCutterUndoRecall
%
% ADR 2008
%
% Part of the new undo mechanism in MClust 3.5.
% Modified from code by JCJ 2007
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

global MClust_Undo MClust_Redo
global MClust_Clusters  MClust_Colors 

if isempty(MClust_Undo)
    warndlg('Nothing to undo.')
    return
end

% push the redo stack
U = MClustCutterUndoGetCurrentState(MClust_Undo{end}.funcname);
MClust_Redo{end+1} = U;

% pop the undo stack
MC_Undo = MClust_Undo{end}; 
MClust_Undo(end) = [];
MClust_Clusters = MC_Undo.clusters;
MClust_Colors = MC_Undo.colors; 

msgbox(['Undid function ', MC_Undo.funcname], 'Undo', 'none', 'non-modal');

global MClust_FilesWrittenYN
MClust_FilesWrittenYN = 'no';

% update tooltips
MClustCutterUndoUpdateTooltip();

% count down autosave
MClustCutterStepAutosave;
