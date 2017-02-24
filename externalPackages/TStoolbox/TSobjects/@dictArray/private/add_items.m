function O = add_items(O, k, v)
% add_item(D, k, v) add items with keys k and values v to dictArray D 
%
% k and v are 
% if key k is already present the corresponding value is replaced
  
% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  ix = key2idx(O, k);
  
  for i = 1:length(ix)
    if ix(i) > 0
      O.values{ix(i)} = v{i};
    else
      O.keys(end+1) = { k{i} };
      O.values(end+1) = { v{i} };
    end
  end
  
  
  
  
  