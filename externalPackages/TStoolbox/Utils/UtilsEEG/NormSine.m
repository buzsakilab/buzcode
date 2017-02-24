function [NS, phase] = NormSine(S)
  
  
  ds = Data(S);
  ts = Range(S, 'ts');
  
  [ns, phase, s1, s2] = NormSineImpl(ds);
  
  NS= tsd(ts(s1:s2), ns');
  phase = tsd(ts(s1:s2), phase');
  