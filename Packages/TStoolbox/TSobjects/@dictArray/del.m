function O = del(O,k)
% D = del(D,k) deletes item with key k in dictArray D
%
% gives an error if k not present

% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  
  if O.keysReadOnly
    error('Keys are not modifiable');
  end

  ix = key2idx(O, k);
  
  if ix == 0
    error('missing key');
  end
  
  keys = O.keys;
  values = O.values;
  l = length(keys);
  keys = keys([1:(ix-1) (ix+1):l]);
  values = values([1:(ix-1) (ix+1):l]);
  
  O.keys = keys;
  O.values = values;