function [D] = Dispositions(N, k);
  
  
  nd = factorial(N) / factorial(N-k);
  
  
  D = zeros(nd, k);
  
  C = gsl_combinations(N, k);
  P = gsl_permutations(k);
  
  nc = size(C, 1);
  np = size(P, 1);
  
  l = 0;
  for i = 1:nc
    for j = 1:np
      l = l+ 1;
      cc = C(i,:);
      D(l,:) = C(i,P(j,:));
    end
  end
  
  