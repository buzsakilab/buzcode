function index = binsearch(data, key)

start = 1;
ends = length(data);
while (start < (ends-1))

      mid = floor((start + end)/2);

      if ((key) == data(mid))
	start = mid;
        ends = mid;
      if ((key) < data(mid))
	ends = mid;
      if ((key) > data(mid))
	start = mid;
end

  if(((key) - data(start)) > (data(start+1) - (key)))
    start = start+1;
  
    index  = start;
    
  