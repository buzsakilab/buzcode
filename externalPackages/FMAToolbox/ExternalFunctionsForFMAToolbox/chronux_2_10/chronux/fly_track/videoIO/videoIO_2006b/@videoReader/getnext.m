function frame = getnext(vr)
%FRAME=GETNEXT(VR)
%  This method is a shortcut for
%    if (next(vr))
%      frame = getframe(vr);
%    else
%      frame = [];
%    end
%
%  Although this method can be more convenient than calling NEXT and GETFRAME
%  separately, there is no way to distinguish between a zero-sized frame and
%  reading past the end of the stream.  
%
%  Typical usage:
%    vr = videoReader('numbers.uncompressed.avi');
%    info = getinfo(vr);
%    for i=1:info.numFrames
%      img = getnext(vr);
%      imshow(img);
%      pause(1/info.fps);
%    end
%    vr = close(vr);
%
%SEE ALSO
%  videoReader
%  videoReader/getframe
%  videoReader/next
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

if (next(vr))
  frame = getframe(vr);
else
  frame = [];
end


