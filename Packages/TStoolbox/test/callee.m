function callee()

  k = evalin('caller', 'aaa;');
  
  display(['in caller k is ' num2str(k)]);
  
  
  c = callingFunction
  
  display(['caller is ' c]);