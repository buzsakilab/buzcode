function [pval, F, df1, df2] = OneWayAnova(X, G)
% [pval, F, df1, df2] = OneWayAnova(X, G)
%
% Computes one-way ANOVA (Between subjects)
% 
% INPUTS:
% X = a vector containing the data points
% G = a vector (same dims as X) containing indices to levels
% (experimental group), 1 corresponds to first level, 2 to second level,
% etc.
%
% OUTPUTS:
% pval: p value
% F: The computed F statistics
% df1: numerator degrees of freedom
% df2: denominator degrees of freedom
  
  
  
  if(any(size(X) ~= size(G)))
    error('X and G must have same size');
  end
  
  
  if any(floor(G) ~= G | G <=0)
    error('G must be a vector of positive integers');
  end
  
  
  
  N = length(G);
  uG = unique((G(:))');
  k = length(uG);
  
  SXT = sum(X);
  SX2T = sum(X.*X);
  
  SST = SX2T - SXT * SXT / N;
  
  SSBG = 0;
  for i = uG
    ix = find(G==i);
    XG = X(ix);
    SSBG = SSBG + (sum(XG))^2 / length(ix);
  end
  
  SSBG = SSBG - SXT * SXT / N;
  
  SSWG = SST - SSBG;
  
  dfBG = k - 1;
  dfWG = N - k;
  
  MSBG = SSBG / dfBG;
  MSWG = SSWG / dfWG;
  
  F = MSBG/MSWG;
  df1 = dfBG;
  df2 = dfWG;
  
  pval = 1 - fp(F, df1, df2);
  
  
  
  
  
  

  