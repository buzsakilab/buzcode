function vr = close(vr)
%VR=CLOSE(VR)
%  Closes video VR and releases any system resources necessary to access it
%  (e.g. threads, file handles, etc.).  Do NOT just clear a videoReader
%  object without first closing its handle:
%
%    % BAD code--typically will leak system resources
%    vr = videoReader(...);
%    ...
%    clear vr; % leaks resources
%
%    % GOOD code
%    vr = videoReader(...);
%    vr = close(vr);
%    clear vr; % okay, but not needed
%
%  After calling CLOSE, VR should not be used any more.
%    vr = videoReader(...);
%    vr = close(vr);
%    next(vr); % BAD
%    vr = videoReader(...);
%    close(vr); % BAD: should reassign result to vr to be safe
%    next(vr); % BAD
%  
%SEE ALSO
%  videoReader
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

feval(vr.plugin, 'close', vr.handle);
vr.handle = nan;
