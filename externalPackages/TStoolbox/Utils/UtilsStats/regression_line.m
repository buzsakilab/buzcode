function [a, b, r, p] = regression_line(X, Y)
%  [a, b, r, p] = regression_line(X, Y)

  ix = find(isfinite(X.*Y));
  
  if(length(ix) < 3)
    a= NaN;
    b = NaN;
    r = NaN;
    p = NaN;
    return
  end
  
  
  
  b = (mean((X(ix)).*(Y(ix))) - mean((X(ix))) * mean((Y(ix)))) / ...
      (mean((X(ix)).^2) - mean((X(ix)))^2);
  
  a = mean((Y(ix))) - b * mean((X(ix)));  
  
  r = corrcoef(X(ix), Y(ix));
  r = r(1,2);
  df = length(ix)-2;
  t_val = r * sqrt(df) / sqrt(1 - r^2);
  p  = 2 * (1 - tcdf(abs(t_val),df));
