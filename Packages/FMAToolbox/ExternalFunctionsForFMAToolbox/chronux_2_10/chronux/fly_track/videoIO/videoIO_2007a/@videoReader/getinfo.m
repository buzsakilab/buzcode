function info = getinfo(vr)
%INFO=GETINFO(VR)
%  Returns a structure whose fields contain information about the opened
%  video object.  The (minimum) set of fields in the INFO structure is
%  shown below:
%    url        String specifying the data source, in the format preferred
%               by the plugin being used.  Sometimes this will be a true
%               URL, sometimes it will be a filename.
%
%    fps        Non-negative number indicating the number of frames per
%               second.  
%
%    height     Integer indicating the height of the video frames in
%               pixels.   
%
%    width      Integer indicating the width of the video frames in pixels.
%
%    numFrames  Integer indicating an estimate of the total number of 
%               frames in the video.  For typical videos, this number is
%               exact.  Users may attempt to read more than numFrames
%               frames at their own risk.  If nHiddenFinalFrames is
%               non-zero, this will typically fail (meaning next/step/seek
%               will return 0) or worse, corrupted data such as an
%               all-black frame may be returned by the codec.  Some plugins
%               and/or their codecs do not supply this information.  If the
%               number of frames is unknown, a negative number is returned.
%               Older ffmpeg versions (notably version 0.4.9-pre1) do not
%               supply this number.
%
%    fourcc     4- or fewer-character string roughly indicating the codec
%               used encode the video.  See http://www.fourcc.org for
%               additional background information and an extensive, but
%               non-comprehensive list of FourCC codes.
%
%    nHiddenFinalFrames
%               Non-negative integer.  Many codecs make it difficult or 
%               impossible to read the last few frames of a file.  When 
%               videoReader thinks that the last few cannot be read, it 
%               automatically guesses how many frames cannot be read,
%               records this number as nHiddenFinalFrames, and sets
%               numFrames to be the number of frames the file claims to
%               contain minus nHiddenFinalFrames.  An individual
%               videoReader plugin (like the ffmpegPopen2 plugin) may choose
%               to allow the user to try reading the frames that might be
%               hidden or it may choose not to allow even trying to read
%               them (like the DirectShow plugin).
%
%  Due to limitations in some file formats, it is not always possible to
%  determine all of these values (or sometimes they are not constant).  In
%  these cases, numerical values are given a value of NaN and string values
%  are blank.
%
%SEE ALSO
%  videoReader
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

%info = feval(vr.plugin, 'getinfo', vr.handle);

[names, vals] = feval(vr.plugin, 'getinfo', vr.handle);
info = cell2struct({vals{:}}, {names{:}}, 2);
