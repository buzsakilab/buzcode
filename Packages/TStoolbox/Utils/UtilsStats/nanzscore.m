function S = nanzscore(X)
  
if min(size(X)) == 1;
    X = X(:);
end

S = NaN(size(X));

for ii=1:size(X,2)
    ix = ~isnan(X(:,ii));
    S(ix,ii) = zscore(X(ix,ii));
end
  
  