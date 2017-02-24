function O = time_units(inp)
  
  if isa(inp, 'units') % copy constructor
    O = inp;
    return
  end
  
  if ~isa(inp, 'char')
    error('argument must be string');
  end
  
  l = get_lookup;
  
  if ~has_key(l, inp)
    error('unrecognized unit');
  end
  
  l.tag = inp 
  
  l = class(l, 'time_units');
  
    
    