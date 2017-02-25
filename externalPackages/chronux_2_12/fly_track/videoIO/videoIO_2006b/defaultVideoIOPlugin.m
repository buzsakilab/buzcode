function pluginName = defaultVideoIOPlugin
%pluginName = defaultVideoIOPlugin
%  Returns the name of the default videoIO (videoReader/videoWriter) plugin.  
%  The videoReader and videoWriter constructors will use the returned plugin if 
%  none is specified.
%
% SEE ALSO:
%   buildVideoIO
%   videoReader
%   videoWriter
%
%Copyright (c) 2007 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

if ispc
  pluginName = 'DirectShow';
else
  pluginName = 'ffmpegPopen2';
end    
