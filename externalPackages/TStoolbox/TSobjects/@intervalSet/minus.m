function s = minus(a, b)

%  Returns the set difference of two intervalSet (overload of the - operator)
%  
%  	USAGE:
%  	s = and(a, b) 

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  s = diff(a, b);
  
  