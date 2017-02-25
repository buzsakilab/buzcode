function frame = getframe(vr)
%FRAME=GETFRAME(VR)
%  Returns the current frame that has been read in from video VR.  Calling
%  GETFRAME does not automatically advance the current frame so the user
%  must first call NEXT, STEP, or SEEK before calling GETFRAME the
%  first time.  This makes it easy to write loops such as:
%    vr = videoReader(...);
%    while (next(vr))
%      img = getframe(vr);
%      ...do something...
%    end
%    vr = close(vr);
%
%SEE ALSO
%  videoReader
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

frame = feval(vr.plugin, 'getframe', vr.handle);

