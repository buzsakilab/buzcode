function p = percentile(V, x);
% p = percentile(V, x);
% 
% computes x-th percentile of vector V
  
  
  
if x > 1 % interpret as percentage
  x = x / 100;
end

  
Vsort = sort(V);
lv = length(V);
ix = floor(lv*x);
if ix < 1
  ix = 1;
end

if x > lv 
  ix = lv;
end



p = Vsort(ix);


