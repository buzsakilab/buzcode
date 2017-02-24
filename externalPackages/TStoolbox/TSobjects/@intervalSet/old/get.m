function val = get(O,prop_name)
% GET Get intervalSet properties from the specified object
% and return the value

% batta 2001
% starting version



switch prop_name
 case 'start'
  val = O.start;
 case 'stop'
  val = O.stop;
 case 'units';
  val = O.units;
end

