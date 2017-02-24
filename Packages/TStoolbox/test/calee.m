function calee()

  k = evalin('caller', 'aaa;');
  
  display(['in caller k is ' num2str(k)]);
  