function MClustCutterUndoStore(U)

% MClustCutterUndoStore(undoStateStruct)
% MClustCutterUndoStore(funcname)
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

global MClust_Undo MClust_Redo MClust_MaxUndoNmbr

if ischar(U)
    U = MClustCutterUndoGetCurrentState(U); % implements funcname case
end
    
% STORE IT!
MClust_Undo{end+1} = U;

% Clear redo stack - it's no longer valid.
MClust_Redo = {};

% Keep length of Undo cue limited to save memroy - JCJ Oct 2007
if (length(MClust_Undo) > MClust_MaxUndoNmbr) 
    MClust_Undo(1) = [];
end

global MClust_FilesWrittenYN
MClust_FilesWrittenYN = 'no';

% update tooltips
MClustCutterUndoUpdateTooltip();

% count down autosave
MClustCutterStepAutosave;
