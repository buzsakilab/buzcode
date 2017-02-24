function display(O)
%   display(O) deals with displaying dictArrays
  k = O.keys;
  v = O.values;
  format compact
  fprintf(1, '\n{');
  
  for i = 1:length(k)
    fprintf(1, '  %s:  ', k{i});
    if ~isempty(v{i})
      disp(v{i});
    else
      fprintf(1, '\n');
    end
    
    fprintf(1, ' ');
  end
  
  fprintf(1, '}\n\n');
  
  
    