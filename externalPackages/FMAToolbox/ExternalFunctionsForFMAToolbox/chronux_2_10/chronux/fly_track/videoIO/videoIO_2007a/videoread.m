function mov=videoread(fname,varargin)
%VIDEOREAD Read an video file in a mannter compatible with Matlab's
%          built-in AVIREAD function, but using a videoReader object to do
%          the work instead. 
%
%   MOV=VIDEOREAD(FNAME)
%     Read the video file FNAME into the Matlab movie structure MOV.  MOV
%     has two fields, "cdata" and "colormap".  Since all videoReader
%     plugins currently output 24-bit images in all cases, colormap will
%     always be blank.    
%
%   MOV=VIDEOREAD(FNAME,VRARGS)
%     The same as MOV=VIDEOREAD(FNAME), but VRARGS are passed to the 
%     videoReader constructor.  VRARGS should be a cell array of arguments.
%     E.g. 
%         mov = videoread('tests/numbers.wmv3.avi', {'DirectShow'});
%
%   MOV=VIDEOREAD(FNAME,INDEX)
%     Reads only the frame(s) specified by INDEX.  INDEX can be a single
%     index or an array of indices into the video stream, where the first
%     frame is  number one.  IMPORTANT NOTES: If you use videoReader
%     directly, it assumes  frames start at 0, not 1.  We use 1-indexing
%     here in VIDEOREAD to maximize  compatibility with MathWork's AVIREAD.  
%     Also, if you find yourself making several calls to VIDEOREAD to 
%     extract different portions of the file and/or if the number of 
%     frames is large, consider using videoReader directly instead of this 
%     wrapper function.
%
%   MOV=VIDEOREAD(FNAME,VRARGS,INDEX)
%     The same as MOV=VIDEOREAD(FNAME,INDEX), but VRARGS are passed to the 
%     videoReader constructor.  VRARGS should be a cell array of arguments.  
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially 
%when using this library on GNU/Linux). 

if nargin==1,
  vrargs = {};
elseif nargin==2,
  if (iscell(varargin{1})),
    vrargs = varargin{1};
  else
    vrargs = {};
    index  = varargin{1};
  end
elseif nargin==3,
  vrargs = varargin{1};
  index  = varargin{2};
end

vr = videoReader(fname, vrargs{:});
if ~(exist('index','var')),
  info = getinfo(vr);
  index = (0:info.numFrames-1); % 0-indexed
else
  index = index - 1; % 0-indexed
end

mov = struct;
currFrame = -1;
if (size(index)>0)
  % index supplied or deduced from file
  for i=1:length(index)
    f = index(i);
    if (f == currFrame+1)
      if ~next(vr), error('Could not advance to frame %d\n', f+1); end
    else
      if ~seek(vr, f), error('Could not seek to frame %d\n', f+1); end
    end
    currFrame = f;
    mov(i).cdata = getframe(vr);
    mov(i).colormap = [];
  end
else
  % no index was supplied and we could not deduce the number of frames 
  % in the video.
  i=1;
  while next(vr)
    mov(i).cdata = getframe(vr);
    mov(i).colormap = [];
    i = i+1;
  end
end

vr = close(vr);

