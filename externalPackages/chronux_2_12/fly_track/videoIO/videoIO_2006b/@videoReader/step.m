function fn = step(vr, delta)
%WORKED=STEP(VR,DELTA)
%  Moves the frame counter by DELTA frames for video VR.  This is a 
%  generalization of NEXT.  Returns 0 on an unsuccessful STEP.  Note that 
%  not all plugins support stepping, especially with negative numbers.  In 
%  the following example, both IM1 and IM2 should be the same for most 
%  plugins.
%    vr = videoReader(...myurl...);
%    if (~next(vr)), error('couldn''t read first frame'); end
%    im1 = getframe(vr);
%    if (~step(vr,-1)), error('could not step back to frame 0'); end
%    im2 = getframe(vr);
%    if (any(im1 ~= im2)), 
%      error('first frame and frame 0 are not the same'); 
%    end
%    vr = close(vr);
%  FNUM should be an integer.
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

fn = feval(vr.plugin, 'step', vr.handle, delta);

