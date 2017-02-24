function [z, p] = binomial_sign_test(i1, n1, i2, n2)
% [z, p] = binomial_sign_test(i1, n1, i2, n2)

  p = (i1+i2)/(n1+n2);
  
  z1 = abs((i1 - n1 *p) / sqrt(n1 * p * (1-p)));
  z2 = abs((i2 - n2 *p) / sqrt(n2 * p * (1-p)));
  
  z = min(z1, z2);
  p = erfc(z/sqrt(2));