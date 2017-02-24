function it = items(O)
% returns cell arry of key/value pairs in O
  
% copyright (c) 2004 Francesco P. Battaglia 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


keys = O.keys;
values = O.values;

  
l = length(keys);

it = cell(l,1);

for i = 1:l
  it{i} = { keys{i}, values{i}  };
end


  
  
  