function [anvA, anvB, anvAB] = MixedTwoWayAnova(X1, X2)
  
  
  
  p = 2; % number of levels of the between subject factor, for this
         % function is limited to 2, represented by the matrices X1, X2
	 
  q = size(X1, 2); % number of levels of the within subjects factor,
                   % represented by the columns of the matrices
		   
  if size(X2, 2) ~= q
    error('X1, X2 must have the same number of columns.');
    
  end
  
  n1 = size(X1, 1); % number fo subjects in X1;
  n2 = size(X2, 1); % number of subjects in X2
  
  N = (n1+n2)*q;
  
  XS = sum(sum(X1.*X1)) + sum(sum(X2.*X2));
  T = ( sum(sum(X1)) + sum(sum(X2)))^2 / N;
  
  A = (sum(sum(X1)))^2 / (n1*q) + sum(sum(X2))^2 / (n2*q);
  
  X12 = [X1; X2];
  
  B = sum(sum(X12).^2 / (n1+n2));
  
  AB = sum(sum(X1).^2)/ n1 + sum(sum(X2).^2) / n2;
  
  AS = sum( sum(X12, 2).^2 ) / q; 

  SA = A - T;
  dfA = p -1; % = 1
  MSA = SA/dfA;
  
  SSWG = AS - A;
  dfSWG = (n1-1) + (n2-1);
  MSSWG = SSWG/dfSWG;
  
  SWS = XS-AS;
  dfWS = (n1+n2)*p*(q-1);
  
  SB = B - T;
  dfB = q-1;
  MSB = SB/dfB;
  
  SAB = AB - A - B + T;
  dfAB = (p-1)*(q-1);
  MSAB = SAB/dfAB;
  
  SBSWG = XS - AB -AS + A;
  dfBSWG = (q-1) * (n1-1+n2-1);
  MSBSWG = SBSWG / dfBSWG;
  
  

  anvA.F = MSA/MSSWG;
  anvA.df1 = dfA;
  anvA.df2 = dfSWG;
  anvA.pval = 1 - fp(anvA.F, anvA.df1, anvA.df2);
  
  
  anvB.F = MSB/MSBSWG;
  anvB.df1 = dfB;
  anvB.df2 = dfBSWG;
  anvB.pval = 1 - fp(anvB.F, anvB.df1, anvB.df2);
  
  anvAB.F = MSAB/MSBSWG;
  anvAB.df1 = dfAB;
  anvAB.df2 = dfBSWG;
  anvAB.pval = 1 - fp(anvAB.F, anvAB.df1, anvAB.df2);
  
  