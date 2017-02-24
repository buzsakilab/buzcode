function vr = videoReader(url, varargin)
% videoReader class constructor
%   Creates a object that reads video streams.  We use a plugin 
%   architecture in the backend to do the actual reading.  For example, 
%   on Windows, DirectShow will typically be used and on Linux, the 
%   ffmpeg library is often used.
%
%   vr = videoReader(url)
%     Opens the given video file for reading using the default plugin.
%     On Windows, 'DirectShow' is used by default and on Linux,
%     'ffmpegPopen2' is used by default.  For most plugins, the url will 
%     really be a filename.
%
%   vr = videoReader(url, ..., 'plugin',pluginName, ...)
%   vr = videoReader(url,pluginName, ...)
%     Opens the file using the manually specified plugin implementation.  
%     Available plugins include:
%
%     'DirectShow': preferred method on Windows
%       - Only available on Windows
%       - See INSTALL.dshow.txt for installation instructions
%       - Can load virtually any video file that can be played in
%         Microsoft's Windows Media Player.  Note that the correct codec
%         must be installed to read a file.  For example, to read 
%         tests/numbers.3ivx.avi, the user must have installed an MPEG4
%         codec such as 3ivx (www.3ivx.com), DivX (www.divx.com), or XviD 
%         (www.xvid.org).
%       - The URL parameter should be a filename.
%         - As a convenience, all forward slashes ('/') are automatically  
%           converted to backslashes ('\')
%
%     'ffmpegPopen2': safe method on Linux
%       - Only supported on GNU/Linux (might work on BSD systems too like Mac 
%         OS X, but this is untested)
%       - See INSTALL.ffmpeg.txt for installation instructions
%       - Creates a separate server process to communicate with the
%         ffmpeg libraries.  
%         - Works when the system's version of GCC is very different from
%           the one that MathWorks used to compile Matlab.
%         - Isolates ffmpeg errors so they typically cannot crash
%           Matlab.  
%         - May allow for more flexible distribution terms for your code 
%           when it uses videoIO (ffmpeg may be compiled with either 
%           the LGPL or GPL license).
%
%     'ffmpegDirect': low-overhead method on Linux
%       - same as ffmpegPopen2, but the ffmpeg libraries are loaded
%         directly by the MEX file.
%         - May not work if MathWorks' and your version of GCC are
%           incompatible. 
%         - Slightly faster than ffmpegPopen2 since there is no
%           interprocess communication overhead.
%
%   vr = videoReader(url, ..., param,arg,...)
%     Allows the user to pass extra configuration arguments to plugin.
%     Currently no plugin arguments are supported right now.  In the 
%     future, we may allow the user to do things like have DirectShow
%     automatically convert to grayscale, or give options to trade off 
%     speed with seeking precision.
% 
% Once you have created a videoReader object, you must next call NEXT,
% SEEK, or STEP at least once so that it will read some frame from disk.  
% *After* calling one of these methods, you may call GETFRAME as many 
% times as you would like and it will decode and return the current frame 
% (without advancing to a different frame).  GETINFO may be called at any 
% time (even before NEXT, SEEK, or STEP).  It returns basic information 
% about the video stream.  Once you are done using the videoReader, make 
% sure you call CLOSE so that any system resources allocated by the plugin 
% may be released.  Here's a simple example of how you might use
% videoReader: 
%
%   % take us to the videoReader directory since we know there's a video
%   % file there.  
%   chdir(fileparts(which('videoReaderWrapper.cpp')));
%
%   % Construct a videoReader object
%   vr = videoReader('tests/numbers.uncompressed.avi');
%
%   % Do some processing on the video and display the results
%   avgIntensity = [];
%   i = 1;
%   figure;
%   while (next(vr))
%     img = getframe(vr);  
%     avgIntensity(i) = mean(img(:));
%     subplot(121); imshow(img);        title('current frame');
%     subplot(122); plot(avgIntensity); title('avg. intensity vs. frame');
%     drawnow;      pause(0.1);         i = i+1;
%   end
%   vr = close(vr);
%
% SEE ALSO:
%   buildVideoMex
%   videoReader/close
%   videoReader/getframe
%   videoReader/getinfo
%   videoReader/getnext
%   videoReader/next
%   videoReader/seek
%   videoReader/step
%   videoWriter
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

if (mod(length(varargin),2) == 0)
  plugin     = defaultVideoIOPlugin;
  pluginArgs = varargin;
else
  plugin     = varargin{1};
  pluginArgs = {varargin{2:end}};
end
[plugin,pluginArgs] = parsePlugin(plugin, pluginArgs);

vr = struct('plugin',mexName(plugin), 'handle',int32(-1));
vr = class(vr, 'videoReader');
[pathstr, name, ext, versn] = fileparts(url);
strArgs = cell(size(pluginArgs));
for i=1:numel(pluginArgs), strArgs{i} = num2str(pluginArgs{i}); end
vr.handle = feval(vr.plugin, 'open', vr.handle, ...
  fullfile(pathstr,[name ext versn]), strArgs{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = mexName(plugin)
n = ['videoReader_' plugin];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plugin,pluginArgs] = parsePlugin(plugin, pluginArgs)
if (length(pluginArgs) > 0)
  [pluginSpecified,idx] = ismember('plugin', {pluginArgs{1:2:end}});
  if pluginSpecified
    plugin = pluginArgs{idx*2};
    pluginArgs = { pluginArgs{1:idx*2-2}, pluginArgs{idx*2+1:end} };
  end
end
