function vw = videoWriter(url, varargin)
% videoWriter class constructor
%   Creates a object that writes video files.  We use a plugin 
%   architecture in the backend to do the actual writing.  For example, 
%   on Windows, DirectShow will typically be used and on Linux, the 
%   ffmpeg library is often used.
%
%   vw = videoWriter(url)
%     Opens the given video file for writing using the default plugin.
%     On Windows, 'DirectShow' is used by default and on Linux,
%     'ffmpegPopen2' is used by default.  For most plugins, the url will 
%     really be a filename.
%
%   vw = videoWriter(url,..., 'plugin',pluginName, ...)
%   vw = videoWriter(url,pluginName)
%     Opens the file using the specified plugin implementation.
%     Available plugins include:
%
%     'DirectShow': preferred method on Windows
%       - Only available on Windows
%       - See INSTALL.dshow.txt for installation instructions
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
%   vw = videoWriter(url, ..., param,arg,...)
%     Allows the user to pass extra configuration arguments to plugin.
%     At present, all parameter names are case sensitive (but in the
%     future they may become case-insensitive).  
%
%     The following parameters are supported by current plugins:
%
%                        Plugin 
%     Parameter      ffmpeg* DShow  Implementation Notes
%     ---------      ------- -----  -----------------------------
%     width             X      X    Width of the encoded video.  Most
%                                   codecs require width to be divisible
%                                   by 2, 4, or 8.  Most users will want
%                                   to explicitly pass this parameter.
%                                   The addframe method will
%                                   automatically resize any images
%                                   according to the value chosen here
%                                   (or a default value if none is
%                                   specified here). 
%
%     height            X      X    Height of the encoded video.  Most
%                                   codecs require height to be divisible
%                                   by 2, 4, or 8.  Most users will want
%                                   to explicitly pass this parameter.
%                                   The addframe method will
%                                   automatically resize any images
%                                   according to the value chosen here
%                                   (or a default value if none is
%                                   specified here). 
%
%     codec             X      X    A string specifying the encoder to
%                                   use.  The exact set of possible
%                                   codecs is highly system-dependent.
%                                   Most users will want to explicitly
%                                   pass this parameter.  To see a list
%                                   of available codecs on a specific
%                                   machine, run:
%                                     codecs = videoWriter([],'codecs');
%
%     fourcc                   X    For the DirectShow plugin, this is a 
%                                   synonym for 'codec'.
%
%     fps               X      X    Frame rate of the recorded video.  
%                                   Note that some codecs only work with 
%                                   some frame rates. 15, 24, 25, 29.97,
%                                   and 30 should work with most codecs.  
%
%     framesPerSecond   X      X    An alias for fps.
%
%     fpsNum, fpsDenom  X      X    This pair of parameters allows frames
%                                   per second to be specified as a
%                                   rational number.  Either both or
%                                   neither parameter must be given.
%
%     framesPerSecond_num           Alias for fpsNum, fpsDenom pair.
%     framesPerSecond_denom
%                       X      X
%
%     bitRateTolerance  X           For supporting codecs, the actual
%                                   bit rate is allowed to vary by +/-
%                                   this value. 
%
%     showCompressionDialog    X    If true (a non-zero number), a dialog
%                                   box is presented to the user allowing
%                                   precise manual selection of the codec
%                                   and its parameters.  Note: sometimes 
%                                   the dialog does not received focus
%                                   automatically so you'll need to 
%                                   ALT-TAB to get to it.
%
%     codecParams              X    A MIME Base64-encoded string describing
%                                   the codec setup parameters for a 
%                                   DirectShow codec.  The contents of this
%                                   string are very codec-specific.  Often,
%                                   The best ways to come up with a string
%                                   like this are to first create a
%                                   videoWriter with the
%                                   'showCompressionDialog' option enabled,
%                                   choose the desired settings, then use
%                                   the GETINFO method to extract the
%                                   'codecParams' value.  Note that this
%                                   representation is the same as used by
%                                   VirtualDub 1.6 and 1.7 in its Sylia
%                                   Script files.  Nearly all useful
%                                   DirectShow codecs can only be
%                                   configured with 'codecParams' and they
%                                   ignore the separate 'bitRate' and
%                                   'gopSize' parameters given below.
%
%     bitRate           X      x    Target bits/sec of the encoded video.
%                                   Supported by most ffmpeg codecs.
%                                   To see whether a particular codec uses
%                                   the bitRate parameter, run the
%                                   testBitRate function in the tests/
%                                   subdirectory (NOTE: very few DirectShow
%                                   codecs support it).   
%
%     gopSize           X      x    Maximum period between keyframes.  GOP
%                                   stands for  "group of pictures" in MPEG
%                                   lingo.  Supported by most ffmpeg
%                                   codecs.  To see whether a particular
%                                   codec uses the gopSize parameter, run
%                                   the testGopSize function in the tests/
%                                   subdirectory (NOTE: very few DirectShow
%                                   codecs support it).
%
%     maxBFrames        X           For MPEG codecs, gives the max
%                                   number of bidirectional frames in a
%                                   group of pictures (GOP).
%
%  codecs = videoWriter([],'codecs')
%  codecs = videoWriter([],pluginName,'codecs')
%  codecs = videoWriter([],'codecs','plugin',pluginName)
%    Queries the backend for a list of the valid codecs that may be used
%    with the 'codec' plugin parameter.  
%
% Once you are done using the videoWriter, make sure you call CLOSE so
% that any system resources allocated by the plugin may be released.
% Here's a simple example of how you might use videoWriter to create
% a video of continually adding more motion blur to an image:
%
%   % Construct a videoWriter object
%   vw = videoWriter('writertest.avi', ...
%                    'width',320, 'height',240, 'codec','xvid');
%   img = imread('peppers.png');
%   h = fspecial('motion',10,5);
%   for i=1:100
%     addframe(vw, img);
%     img = imfilter(img, h);
%   end
%   vw=close(vw);
%
% SEE ALSO:
%   buildVideoMex
%   videoWriter/addframe
%   videoWriter/close
%   videoReader
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

if (numel(url)==0)
  % static method call
  if (mod(length(varargin),2) == 0)
    plugin       = varargin{1};
    staticMethod = varargin{2};
    methodArgs   = {varargin{3:end}};
  else
    plugin       = defaultVideoIOPlugin;
    staticMethod = varargin{1};
    methodArgs   = {varargin{2:end}};
  end
  [plugin,methodArgs] = parsePlugin(plugin, methodArgs);

  vw = feval(mexName(plugin), staticMethod, int32(-1), methodArgs{:});
  
else
  % constructor call
  if (mod(length(varargin),2) == 0)
    plugin     = defaultVideoIOPlugin;
    pluginArgs = varargin;
  else
    plugin     = varargin{1};
    pluginArgs = {varargin{2:end}};
  end
  plugin = parsePlugin(plugin, pluginArgs);
  
  vw = struct('plugin',mexName(plugin), 'handle',int32(-1), ...
              'w',int32(-1), 'h',int32(-1));
  vw = class(vw, 'videoWriter');
  [pathstr, name, ext, versn] = fileparts(url);
  strArgs = cell(size(pluginArgs));
  for i=1:numel(pluginArgs), strArgs{i} = num2str(pluginArgs{i}); end
  [vw.handle,vw.w,vw.h] = feval(vw.plugin, 'open', vw.handle, ...
                                fullfile(pathstr,[name ext versn]), ...
                                strArgs{:});
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = mexName(plugin)
n = ['videoWriter_' plugin];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plugin,pluginArgs] = parsePlugin(plugin, pluginArgs)
if (length(pluginArgs) > 0)
  [pluginSpecified,idx] = ismember('plugin', {pluginArgs{1:2:end}});
  if pluginSpecified
    plugin = pluginArgs{idx*2};
    pluginArgs = { pluginArgs{1:idx*2-2}, pluginArgs{idx*2+1:end} };
  end
end
