function S = sem(X)
  
if min(size(X)) == 1;
    X = X(:);
end

S = std(X) / sqrt(size(X,1));
  
  