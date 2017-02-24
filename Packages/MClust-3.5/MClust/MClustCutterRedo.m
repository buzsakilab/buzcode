function GeneralizedCutterRedo()
% GeneralizedCutterRedo - JCJ Oct 2007

global MClust_Undo MClust_Clusters  MClust_Colors 
global MClust_UndoPtr %JCJ Oct 2007


% Increment Undo pointer if further undo information exists and last non-undid state exists.
if ~exist('MClust_UndoPtr','var') || isempty(MClust_UndoPtr) || ~(MClust_UndoPtr<(length(MClust_Undo)-1))
    warndlg('No steps to redo.')
    return
else
    MClust_UndoPtr=MClust_UndoPtr+1;
end

% Use state before current undo (i.e. undo+one step).
MC_Redo = MClust_Undo(MClust_UndoPtr+1); % MClust_UndoPtr index ~JCJ Oct 2007
% GeneralizedCutterStoreUndo('Undo'); ~ commented out by JCJ Oct 2007
MClust_Clusters = MC_Redo.clusters;
MClust_Colors = MC_Redo.colors; % added restore of color information - JCJ Sept 2007

% function name redone is store under current undo
msgbox(['Re-applied function ', MClust_Undo(MClust_UndoPtr).funcname], 'Redo', 'none', 'modal');


clear MC_Undo
