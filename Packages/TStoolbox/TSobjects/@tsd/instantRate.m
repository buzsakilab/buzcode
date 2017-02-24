function ir = instantRate(tsa, TimeUnits)

%  Computes instantaneous rate as the inverse of inter-event interval
%  	
%  	USAGE:
%  	ir = instantRate(tsa, TimeUnits)
%  	
%  	INPUTS:
%  	tsa 	  - a tsd object
%  	TimeUnits - a units object or the abbreviation string
%  	
%  	OUTPUT:
%  	ir - a tsd object where the data contains the instantaneous rate in the
%  	     units specified by TimeUnits (to the -1)

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  

  
  if nargin == 1
    TimeUnits = time_units('s');
  else 
    TimeUnits = units(TimeUnits);
  end
  
  
  t = Range(tsa, TimeUnits);
  
  ir = tsd(t(1:end-1), 1./diff(t));
  
  