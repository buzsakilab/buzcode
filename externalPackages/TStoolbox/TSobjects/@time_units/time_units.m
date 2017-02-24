function O = time_units(inp)
  
  if isa(inp, 'time_units') % copy constructor
    O = inp;
    return
  end

  if isa(inp, 'units')
    if strcmp(quantity(inp), 'time')
      P = inp;
    else
      error('must be a time unit');
    end
  else
    
    P = units('time', inp);
  end
  
  O.t = 'time_units';  
  O = class(O, 'time_units', P);

  