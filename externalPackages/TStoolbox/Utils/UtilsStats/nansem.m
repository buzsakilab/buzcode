function S = sem(X)
  
if min(size(X)) == 1;
    X = X(:);
end

S = nanstd(X) ./ sqrt(sum(~isnan(X)));
  
  