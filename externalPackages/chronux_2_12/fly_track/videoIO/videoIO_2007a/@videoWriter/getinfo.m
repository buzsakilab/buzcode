function info = getinfo(vw)
%INFO=GETINFO(VW)
%  Returns a structure whose fields contain information about the opened
%  video object.  This structure may be flattened to a cell array and used 
%  with the videoWriter constructor to recreate the current video file.
%
%SEE ALSO
%  videoWriter
%
%Copyright (c) 2006 Gerald Dalley
%See "MIT.txt" in the installation directory for licensing details (especially
%when using this library on GNU/Linux). 

[names, vals] = feval(vw.plugin, 'getinfo', vw.handle);
info = cell2struct({vals{:}}, {names{:}}, 2);

% be nice and convert anything that looks like a number into one.
for i=1:length(names)
  name = names{i};
  if (length(info.(name)) < 64)
    num  = str2num(info.(name));
    if numel(num) > 0
      info.(name) = num;
    end
  end
end
