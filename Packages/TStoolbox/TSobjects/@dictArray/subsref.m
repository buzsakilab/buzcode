function v = subsref(O, S)
% v = subsref(O, S) deals with { } references
% 
% B = O{k}  returns the value corresponding to key k, error if not
% existing key
  
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
  
  if ix == 0
    error('Missing key');
  end
  
  v = O.values{ix};
  