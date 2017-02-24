function worked = next(vr)
%WORKED=NEXT(VR)
%  Advances to the next frame in the video.  Returns 0 if the next frame is
%  not available.  
%
%  After the videoReader constructor is called, NEXT, SEEK, or STEP should
%  be called at least once before GETFRAME is called. 
%
%  See videoReader/getinfo for more details on how many frames can be read
%  from a video. 
%  
%SEE ALSO
%  videoReader
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

worked = feval(vr.plugin, 'next', vr.handle);

