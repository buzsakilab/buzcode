function Cell2List(carray, list)
  
  
  fid = fopen(list, 'w');
  for i = 1:length(carray)
    fprintf(fid, '%s\n', carray{i});
  end
  
  fclose(fid);
  