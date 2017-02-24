function ix = key2idx(O, k);
  
% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
% this waits to be rewritten in C  
  
  if ~isa(k, 'cell')
    k = { k };
  end
  
      
  keys = O.keys;
  
  ix = zeros(length(k), 1);
  
  for i = 1:length(keys)
    for j = 1:length(k)
      if strcmp(k, keys{i})
	ix(j) = i;
	break
      end
    end
  end
  
  
  