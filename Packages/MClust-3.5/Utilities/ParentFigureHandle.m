function H = ParentFigureHandle(objHandle)

% H = ParentFigureHandle(objHandle)
%
% INPUTS
%    objHandle - valid object handle
% 
% OUTPUTS
%    H - figure handle containing object or root if no figure found
%
% ADR 1998
% version L4.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

H = objHandle;
while ~strcmp(get(H, 'Type'),'figure') & ~strcmp(get(H, 'Type'),'root')
   H = get(H, 'Parent');
end