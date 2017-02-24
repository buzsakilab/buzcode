function [pval, F, df1, df2] = RepMeasAnova(X) 
  
  n = size(X, 1);
  k = size(X, 2);
  N = n * k;
  
  XT = sum(X, 1);
  SXT = sum(XT);
  X2T = sum(X.*X, 1);
  SX2T = sum(X2T);
  SXj = sum(X, 2);
  
  
  SST = SX2T - SXT * SXT / N;
  SSBC = sum(XT.*XT) / n - SXT * SXT / N;
  SSBS = sum(SXj.*SXj) / k -  SXT * SXT / N;
  SSres = SST - SSBC - SSBS;
  
  dfBC = k-1;
  dfBS = n-1;
  dfres = (n-1)*(k-1);
  
  MSBC = SSBC / dfBC;
  MSBS = SSBS / dfBS;
  MSres = SSres/dfres;

  F = MSBC / MSres;
  df1 = dfBC;
  df2 = dfres;
  
  pval = 1 - fp(F, df1, df2);
  
  
  
  