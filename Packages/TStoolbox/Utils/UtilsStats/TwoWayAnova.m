function [anvA, anvB, anvAB] = TwoWayAnova(X, GA, GB)
% [pval, F, df1, df2] = OneWayAnova(X, G)
%
% Computes one-way ANOVA (Between subjects)
% 
% INPUTS:
% X = a vector containing the data points
% GA, GB = two vectors (same dims as X) containing indices to levels
% (experimental group), 1 corresponds to first level, 2 to second level,
% etc. for factor A, and factor B
%
% OUTPUTS:
% anvA, anvB, anvAB, strucutures containing the results for main effect
% of A and B, and interaction, with members:  
% pval: p value
% F: The computed F statistics
% df1: numerator degrees of freedom
% df2: denominator degrees of freedom
  
  
  
  if(any(size(X) ~= size(GA) | size(X) ~= size(GB)))
    error('X, GA, GB must have same size');
  end
  
  
  if any(floor(GA) ~= GA | GA <=0)
    error('GA must be a vector of positive integers');
  end
  
   if any(floor(GB) ~= GB | GB <=0)
    error('GB must be a vector of positive integers');
  end
  
 
  
  N = length(GA);
  uGA = unique(GA);
  kA = length(uGA);
  uGB = unique(GB);
  kB = length(uGB);
  
  GG = GA * (kB+1) + GB;

  GH = hist(GG, 1:max(GG));
  
  GH = GH(find(GH));
  
  if any(diff(GH))
    not_equal_size = 1;
  else
    not_equal_size = 0;
  end
  

  nh = kA*kB / sum(1 ./ GH);
  
  if ~not_equal_size
  
  
    SXT = sum(X);
    SX2T = sum(X.*X);

    XS = sum(X.*X);
  
    SST = SX2T - SXT * SXT / N;
    T = (sum(X))^2 / N;
    A = 0;
    for i = uGA
      ix = find(GA==i);
      XG = X(ix);
      A = A + (sum(XG))^2 / length(ix);
    end
  
    B = 0;
    for j = uGB
      ix = find(GB==j);
      XG = X(ix);
      B = B + (sum(XG))^2 / length(ix);
    end
  
    AB = 0;
    for i = uGA
      for j = uGB
	ix = find(GA == i & GB == j);
	if isempty(ix)
	  anvA.pval = NaN;
	  anvA.F = NaN;
	  anvA.df1 = NaN;
	  anvA.df2 = NaN;
	  anvB.pval = NaN;
	  anvB.F = NaN;
	  anvB.df1 = NaN;
	  anvB.df2 = NaN;
	  anvAB.pval = NaN;
	  anvAB.F = NaN;
	  anvAB.df1 = NaN;
	  anvAB.df2 = NaN;
	  return
	end
	XG = X(ix);
	ab(i,j) = (sum(XG))^2 / length(ix);
	AB = AB + ab(i,j);
      end
    end
    SSWG = XS - AB;
  else
    XS = sum(X.*X);
    AB = 0;
    AB_std = 0;
    for i = uGA
      for j = uGB
	ix = find(GA == i & GB == j);
	if isempty(ix)
	  anvA.pval = NaN;
	  anvA.F = NaN;
	  anvA.df1 = NaN;
	  anvA.df2 = NaN;
	  anvB.pval = NaN;
	  anvB.F = NaN;
	  anvB.df1 = NaN;
	  anvB.df2 = NaN;
	  anvAB.pval = NaN;
	  anvAB.F = NaN;
	  anvAB.df1 = NaN;
	  anvAB.df2 = NaN;
	  return
	end
	XG = X(ix);
	ab(i,j) = sum(XG)*nh/length(ix); 
	
	
	AB = AB + (ab(i,j))^2/nh;
	AB_std = AB_std + (sum(XG))^2 / length(ix);
      end
    end
    
    T = (sum(sum(ab)))^2 / (nh * kA * kB);
    
    sa = sum(ab, 2);
    A = sum(sa.^2 / (nh *kB));
    
    sb = sum(ab,1);
    B = sum(sb.^2 / (nh*kA));
    SSWG = XS - AB_std;
  end
  
    
  
    
  
  SST = XS - T;
  SSBG = AB - T;
  SSA = A - T;
  SSB = B - T;
  SSAB = AB - A - B + T;

  
  dfA = kA - 1;
  dfB = kB - 1;
  dfAB = dfA * dfB;
  dfBG = kA*kB - 1;
  dfT = N - 1;
  
  dfWG = dfT - dfBG;
  
  MSA = SSA / dfA;
  MSB = SSB / dfB;
  MSAB = SSAB / dfAB;
  MSWG = SSWG / dfWG;
  


  
  anvA.F = MSA/MSWG;
  anvA.df1 = dfA;
  anvA.df2 = dfWG;
  anvA.pval = 1 - fp(anvA.F, anvA.df1, anvA.df2);
  
  
  anvB.F = MSB/MSWG;
  anvB.df1 = dfB;
  anvB.df2 = dfWG;
  anvB.pval = 1 - fp(anvB.F, anvB.df1, anvB.df2);
  

  anvAB.F = MSAB/MSWG;
  anvAB.df1 = dfAB;
  anvAB.df2 = dfWG;
  anvAB.pval = 1 - fp(anvAB.F, anvAB.df1, anvAB.df2);
  


  

  