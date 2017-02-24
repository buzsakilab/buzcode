function MClustCutterUndoUpdateTooltip()

% MClustUpdateUndoTooltip
%
% ADR 2008
%
% Part of the new undo mechanism in MClust 3.5.
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

global MClust_Undo MClust_Redo

undoObject = findobj('Tag', 'Undo');
redoObject = findobj('Tag', 'Redo');

if isempty(MClust_Undo)
    set(undoObject, 'Tooltip', 'Undo one step [nothing to undo]');
else    
    set(undoObject, 'Tooltip', ...
        sprintf('Undo one step [Step backward]: %s', MClust_Undo{end}.funcname));
end
if isempty(MClust_Redo)
    set(redoObject, 'Tooltip', 'Redo one step [nothing to redo]');
else
    set(redoObject, 'Tooltip', ...
        sprintf('Redo one step [Step forward]: %s', MClust_Redo{end}.funcname));
end