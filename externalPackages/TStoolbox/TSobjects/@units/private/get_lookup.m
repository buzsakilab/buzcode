function l = get_lookup()
  
  ltime = dictArray({ { 'ts', 1},
		  { 'ms', 0.001},
		  { 's', 1},
		  { 'min', 60},
		  { 'h', 3600} });

  ptime = dictArray({ { 'none', 1 } });
    
  l = dictArray({ { 'time', ltime},
		  { 'none', ptime} });
  
  