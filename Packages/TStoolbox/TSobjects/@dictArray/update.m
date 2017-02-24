function O = update(O, A)

% D = update(D, A) updates all the items in D with corresponding items in A  
% 
% items non existing in D will be created anew.  
  
% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


  if O.keysReadOnly
    error('Keys are not modifiable');
  end

  
  k = keys(A);
  v = values(A);
  
  O = add_items(O, k, v);
  
  
  