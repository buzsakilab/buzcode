function hk = has_key(O, k)
% hk = has_key(D, k) returns 1 if key k is present in D, 0 otherwise
  
% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  ix = key2idx(O, k);
  
  hk = ix > 0;
  
  