function O = subsasgn(O, S, B)
% v = subsref(O, S) deals with { } assignments
% 
% O{k} = B  assigns the value B to key k. If key not existing item will
% be created 
  
% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  
  if ~strcmp(S.type, '{}')
    error('reference with { }');
  end
  
  k = S.subs{1};   
  
  if ~isa(k, 'char')
    error('k must be a string');
  end
  
  ix = key2idx(O, k);
  
  if ix > 0
    O.values{ix} = B;
  else
  if O.keysReadOnly
    error('Keys are not modifiable');
  end
    O.keys(end+1) = { k };
    O.values(end+1) = { B } ;
  end
  