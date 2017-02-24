function [z, p] = binomial_test(i, N, p)
  
  z = (i - N *p) / sqrt(N * p * (1-p)); 
  

  p = erfc(z/sqrt(2));