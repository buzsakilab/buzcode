function vw = close(vw)
%VW=CLOSE(VW)
%  Closes video VW and releases any system resources necessary to access it
%  (e.g. threads, file handles, etc.).  Do NOT just clear a videoWriter
%  object without first closing its handle:
%
%    % BAD code--typically will leak system resources
%    vw = videoWriter(...);
%    ...
%    clear vw; % leaks resources, probably results in a corrupted file
%
%    % GOOD code
%    vw = videoWriter(...);
%    ...
%    vw = close(vw);
%    clear vw; % okay, but not needed
%
%  After calling CLOSE, VW should not be used any more.
%    vw = videoWriter(...);
%    vw = close(vr);
%    next(vw); % BAD
%    vw = videoWriter(...);
%    close(vw); % BAD: should reassign result to vw to be safe
%    next(vw); % BAD
%  
%SEE ALSO
%  videoWriter
%
%Copyright (c) 2007 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

feval(vw.plugin, 'close', vw.handle);
vw.handle = nan;
