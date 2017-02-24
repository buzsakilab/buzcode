function clear(O)
% remove all items from dictArray O

% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  
  if O.keysReadOnly
    error('Keys are not modifiable');
  end
  
  O.keys = {};
  O.values = {};
  