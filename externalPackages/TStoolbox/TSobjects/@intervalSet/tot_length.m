function L = tot_length(O, TimeUnits);

%  Returns total length of intervalSet
%  
%  	USAGE:
%  	L = tot_length(O, TimeUnits), 
%  	
%  	Returns the total length in the units specified by TimeUnits
%  	(defaults to ts)

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  if nargin == 1
    TimeUnits = time_units('ts');
  else 
    TimeUnits = time_units(TimeUnits);
  end
  
  
  

  L = sum(End(O, TimeUnits) - Start(O, TimeUnits));
