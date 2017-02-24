function MClustCutterRecallUndo()

% MClustCutterRecallUndo
% added restore of color information - JCJ Sept 2007
global MClust_Undo MClust_Clusters  MClust_Colors 
global MClust_UndoPtr %JCJ Oct 2007


% Decrement Undo pointer - JCJ Oct 2007
if ~exist('MClust_UndoPtr','var') || isempty(MClust_UndoPtr) || MClust_UndoPtr<1
    warndlg('No undo has been stored.')
    return
end

% if stepping back for first time since an action, save current state ~ JCJ Oct 2007
if MClust_UndoPtr==length(MClust_Undo) 
    OldPtr=MClust_UndoPtr;
    MClustCutterStoreUndo(['Undo:' MClust_Undo(MClust_UndoPtr).funcname]); 
    MClust_UndoPtr=OldPtr; 
end

MC_Undo = MClust_Undo(MClust_UndoPtr); % added MClust_UndoPtr index ~JCJ Oct 2007
MClust_Clusters = MC_Undo.clusters;
MClust_Colors = MC_Undo.colors; % added restore of color information - JCJ Sept 2007
msgbox(['Undid function ', MC_Undo.funcname], 'Undo', 'none', 'modal');

MClust_UndoPtr=MClust_UndoPtr-1; %JCJ Oct 2007

clear MC_Undo
