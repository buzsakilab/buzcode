function SE = StartEnd(O, TimeUnits)

%  Starting and Ending times of for each interval in the set
%  
%  	USAGE:
%  	S = StartEnd(O, TimeUnits)
%  	
%  	INPUTS: 
%  	O - an intervalSet object
%  	TimeUnits (optionnal) - a units object or the abbreviation string
%  	
%  	OUTPUT:
%  	SE = an nx2 array of starting and ending points for the n intervals in
%  	the set.  Column 1 are start points, Column2 are end points


% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


  if nargin < 1 | nargin > 2
    error('Call with one or two arguments');
  end
  
  if nargin == 1
    TimeUnits = time_units('ts');
  end  
  


S = Start(O,TimeUnits);
E = End(O,TimeUnits);

SE = cat(2,S,E);