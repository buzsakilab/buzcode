function c = callingFunction
  
  [st, i] = dbstack;
  
  if length(st) < (i+2)
    c = 'base';
  else
    c = st(i+2).name;
    [p, c, e] = fileparts(c);
  end
  