function S = End(O, TimeUnits)
% S = End(O, TimeUnits)%
% INPUTS: 
% O: a ts object
% TimeUnits: a units object or the abbreviation string

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


  if nargin < 1 | nargin > 2
    error('Call with one or two arguments');
  end
  
 if nargin == 1
    TimeUnits = time_units('ts');
  end  
  
  S = O.stop;
  
  if isa(TimeUnits, 'char')
    TimeUnits = time_units(TimeUnits);
  end
  
  cnvrt = convert(O.units, TimeUnits);
  
  if cnvrt ~= 1
    S = S * cnvrt;
  end
  