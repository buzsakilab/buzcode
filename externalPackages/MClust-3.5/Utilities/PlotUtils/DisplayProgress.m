function DisplayProgress(iCounter, maxCounter, varargin)

% DisplayProgress(iCounter, maxCounter, varargin)
% DisplayProgress close
%
% INPUTS
%     iCounter = progress so far
%     maxCounter = when to stop
% 
% PARAMETERS
%   UseGraphics (default true): if true shows bar on screen else prints to stderr
%   Title (default 'Progress so far'): title for waitbar if used
%   EOL (default 80): end of line
%
% ADR 1998
% version L5.2
% status PROMOTED

% v5.0 19 nov 98 Changed parameter from GraphicsTitle to Title
% v5.1 21 jan 99 if only one to counter then skip
% v5.2 31 jan 99 now allows 'close' input
% v5.4 19 Dec 04 now allows close after manually closed wait bar (JCJ)
% v5.5 03 Jan 05 now allows close after manually closed wait bar (JCJ)

% close handle
if ischar(iCounter) & strcmp(iCounter, 'close') 
   global DisplayProgressHandle
   try  % allows for if the wait bar was closed manually -- JCJ Dec 2004
       close(DisplayProgressHandle);
   end
   clear global DisplayProgressHandle
   return
end

if maxCounter == 1; return; end

%-------------------
% parameters
adrlib;

UseGraphics = true;
Title = 'Progress so far';

EOL = 80;
SoFar = 10;

extract_varargin;

%--------------------
if UseGraphics
   global DisplayProgressHandle
   if isempty(DisplayProgressHandle)
      DisplayProgressHandle = waitbar(0, Title);
   else
      waitbar(iCounter/maxCounter);
   end
   if iCounter == maxCounter
       try  % allows for if the wait bar was closed manually -- JCJ Jan 2005
           close(DisplayProgressHandle);
       end
      clear global DisplayProgressHandle
   end
   drawnow;
else
   if iCounter == 1
      fprintf(2, [Title ': .']);  
   elseif iCounter == maxCounter
      fprintf(2, '\n');
   elseif rem(iCounter,EOL) == 0
      fprintf(2, '\n');
   elseif rem(iCounter,10) == 0
      fprintf(2, '!');
   else
      fprintf(2, '.');
   end
end

