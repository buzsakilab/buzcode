function A = subsasgn(A, S, B)

% overload of operator {} for element. 
% if subscripts are numieric, the correspondign elements of the array are returned, if
% subscript is a string, all the tsd's are returned that contain that
% string in their name 
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  if length(S) > 1
    error('subsref defined only for one level of subscripting');
  end
  
  if ~strcmp(S.type, '{}')
    error('only {} subscripting');
  end
  
  
  num_subs = 1;
  for i = 1:length(S.subs)
    if ~isa(S.subs{i}, 'numeric')
      num_subs = 0;
    end
  end
  if num_subs
    A.C(S.subs{:}) = {B};
    return
  end
  
  error('subscript of unrecognized type');
    
      
      
  
  
  
    