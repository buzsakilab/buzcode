function [t, p] = regression2_test(X1, Y1, X2, Y2)
%   [a, b, r, p] = regression_line(X, Y)

  ix1 = find(isfinite(X1.*Y1));
  n1 = length(ix1);
  
  b1 = (mean((X1(ix1)).*(Y1(ix1))) - mean((X1(ix1))) * mean((Y1(ix1)))) / ...
      (mean((X1(ix1)).^2) - mean((X1(ix1)))^2);
  
  a1 = mean((Y1(ix1))) - b1 * mean((X1(ix1)));  
 
  r = corrcoef(X1(ix1), Y1(ix1));
  r1 = r(1,2);

  sX1 = std(X1(ix1));
  sY1 = std(Y1(ix1));
  sYX1 = sY1 * sqrt((1 - r1^2) * (n1-1)/(n1-2));
  
  
  
  ix2 = find(isfinite(X2.*Y2));
  n2 = length(ix2);
  
  b2 = (mean((X2(ix2)).*(Y2(ix2))) - mean((X2(ix2))) * mean((Y2(ix2)))) / ...
      (mean((X2(ix2)).^2) - mean((X2(ix2)))^2);
  
  a2 = mean((Y2(ix2))) - b2 * mean((X2(ix2)));  
 
  r = corrcoef(X2(ix2), Y2(ix2));
  r2 = r(1,2);

  sX2 = std(X2(ix2));
  sY2 = std(Y2(ix2));
  sYX2 = sY2 * sqrt((1 - r2^2) * (n2-1)/(n2-2));
  
  
  
  
  t = (b2-b1) / sqrt((sYX1^2 / (sX1^2 * (n1 - 1))) + ...
		     (sYX2^2 / (sX2^2 * (n2 - 1))) );
  
  
  df = n1+n2-4;
  p = 2 * (1 - tp(abs(t),df));
 