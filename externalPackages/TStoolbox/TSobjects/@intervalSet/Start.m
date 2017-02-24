function S = Start(O, TimeUnits)

%  Starting times of for each interval in the set
%  
%  	USAGE:
%  	S = Start(O, TimeUnits)
%  	
%  	INPUTS: 
%  	O - an intervalSet object
%  	TimeUnits (optionnal) - a units object or the abbreviation string
%  	
%  	OUTPUT:
%  	S = an array of starting points for each interval in the set

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  if nargin < 1 | nargin > 2
    error('Call with one or two arguments');
  end
  
  if nargin == 1
    TimeUnits = time_units('ts');
  end  
  
  S = O.start;
  
  if isa(TimeUnits, 'char')
    TimeUnits = time_units(TimeUnits);
  end
  
  cnvrt = convert(O.units, TimeUnits);
  
  if cnvrt ~= 1
    S = S * cnvrt;
  end
  