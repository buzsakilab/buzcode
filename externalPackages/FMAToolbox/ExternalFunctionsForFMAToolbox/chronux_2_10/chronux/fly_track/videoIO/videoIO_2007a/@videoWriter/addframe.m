function addframe(vw,img)
%WORKED=ADDFRAME(VW)
%  Writes a video frame.  VW must be a videoWriter object and IMG a 2D
%  numeric array with 1 or 3 channels that is autoconverted to a uint8
%  image as follows: 
%
%    type    assumed range
%    ----    -------------
%    uint8   0 to 255
%    double  0 to 1
%    logical 0 to 1
%
%SEE ALSO
%  videoWriter 
%
%Copyright (c) 2007 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

[h,w,d] = size(img);

if (isa(img, 'uint8'))
  if (h ~= vw.h || w ~= vw.w)
    img = uint8(255*imresize(double(img)/255, [vw.h vw.w]));
  else
    % no changes needed
  end
  
elseif (isa(img, 'double') || islogical(img))
  if (h ~= vw.h || w ~= vw.w)
    img = uint8(255*imresize(img, [vw.h vw.w]));
  else
    img = uint8(255*img);
  end
  
else
  error('Invalid image type.');
end

if (d == 1)
  img = repmat(img, [1 1 3]);
end

feval(vw.plugin, 'addframe', vw.handle, img);
