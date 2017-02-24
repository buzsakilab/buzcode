function MClustCutterStoreUndo(funcname)

% MClustCutterStoreUndo(funcname)
% added restore of color information - JCJ Sept 2007
global MClust_Undo MClust_Clusters  MClust_Colors 
global MClust_UndoPtr % JCJ Oct 2007
global MClust_MaxUndoNmbr
% global MClust_ClusterFileNames - no longer used


% Increment Undo pointer - JCJ Oct 2007
if exist('MClust_UndoPtr','var') && ~isempty(MClust_UndoPtr)
    MClust_UndoPtr=MClust_UndoPtr+1;
else
    MClust_UndoPtr=1;
end

MClust_Undo(MClust_UndoPtr).colors = MClust_Colors; % added save of color information - JCJ Sept 2007
MClust_Undo(MClust_UndoPtr).clusters = MClust_Clusters;
MClust_Undo(MClust_UndoPtr).funcname = funcname;
% MClust_Undo(MClust_UndoPtr).clusternames = MClust_ClusterFileNames; % ADR no longer used

%if there are steps ahead of the new Undo step remove them -JCJ Oct 2007
if length(MClust_Undo)>MClust_UndoPtr 
    MClust_Undo((MClust_UndoPtr+1):end)=[];
end

% Keep length of Undo cue limited to save memroy - JCJ Oct 2007
if (MClust_UndoPtr>MClust_MaxUndoNmbr) 
    MClust_Undo(1)=[];
    MClust_UndoPtr=MClust_UndoPtr-1;
end


global MClust_FilesWrittenYN
MClust_FilesWrittenYN = 'no';

% count down autosave
MClustCutterStepAutosave;
