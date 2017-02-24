function DIRSTACK = pushdir(newdir);

% dirstack = pushdir(newdir)
%
% Pushes the current dir onto the directory stack.
% Then cd's to the newdir if given.
% Maintains the directory stack in a global variable.
%
% ADR 1998
% version L4.1
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

global DIRSTACK

if isempty(DIRSTACK)
   DIRSTACK = {pwd};
else
   DIRSTACK = [{pwd} DIRSTACK];
end

if nargin == 1
   cd(newdir);
end
