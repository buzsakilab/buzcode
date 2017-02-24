function fixResources(matfile)
  
  load(matfile, 'resources')
  
  k = keys(resources);
  for i = 1:length(resources)
    R = resources{k{i}};
    R.createdBy = 'Unknown';
    R = Resource(R);
    resources{k{i}} = R;
  end
  
  save(matfile, 'resources', '-append');
  