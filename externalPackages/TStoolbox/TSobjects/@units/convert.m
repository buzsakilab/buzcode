function c = convert(A, B)
% c = convert(A, B) returns multiplicative factor to convert unit A into unit  B
%
% the conversion is such that if a is data expressed in unit A
% b = c * a is the same data expressed in unit B

% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  if(~strcmp(A.quantity, B.quantity))
    error('cannot convert inhomogeneous units');
  end
  
  c = A.value / B.value;
  
  
