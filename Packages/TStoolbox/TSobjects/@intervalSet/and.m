function s = and(a, b)

%  Intersection between two intervalSet objects (overload of the & operator)
%  
%  	USAGE:
%  	s = and(a, b) 

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  s = intersect(a, b);
  
  