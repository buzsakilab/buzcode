function C = ReactCC(X, Y)
% C = ReactCC(X, Y)
%
% coompute the correlation coefficient between matching columns in the
% matrices X and Y, taking care of NaNs
  
  
  if any(size(X) ~= size(Y))
    error('X and Y must have same size');
  end
  mk = zeros(size(X));
  mk(find(isnan(X .* Y))) = NaN;

  X = X + mk;
  Y = Y + mk;
  
  clear mk
  
  mu_x = nanmean(X);
  mu_y = nanmean(Y);
  
  sigma_x = nanstd(X);
  sigma_y = nanstd(Y);
  

  norm_xy = sum(~isnan(X));
  
  X(isnan(X)) = 0;
  Y(isnan(Y)) = 0;
  
  s_xy = sum(X.*Y) ./ norm_xy;
  
  C = (s_xy - (mu_x .* mu_y)) ./ (sigma_x .* sigma_y);
  
  
  
