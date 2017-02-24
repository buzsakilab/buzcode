function [pval, F, df1, df2] = AnovaWithinSubjects(X)



[n,k] = size(X);
N = n*k;
% N the number of subjects  K the number of levels
SS = sum(X,2);
SX = sum(X, 1);
SX2 = sum(X.*X);
SXT = sum(SX);
SX2T = sum(SX2);
SST = SX2T - SXT^2/N;
SSBC = sum(SX.^2)/n - SXT^2/N;
SSBS = sum(sum(SS.^2))/k- SXT^2/N;
SSres = SST-SSBC-SSBS;
dfBC = k-1;
dfBS = n-1;
dfres = (n-1)*(k-1);
MSBC = SSBC/dfBC;
MSres = SSres/dfres;

F = MSBC/MSres;

df1 = dfBC;
df2 = dfres;
  pval = 1 - fp(F, df1, df2);
