function v = get(O, k, varargin)
% v = get(D, k, x) returns item in dictArray D with key k, if key present, otherwise x 
%
% it doesn't give an error if key is not present. In that case, if x is
% not passed, it returns the empty array  
  
% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
x = [];

if length(varargin) > 1
  error('call as v = get(D, k[, x])')
end

if length(varargin) > 0
  x = varargin{1};
end

ix = key2idx(O, k);

if ix == 0
  v = x;
else
  v = O.values{ix};
end


