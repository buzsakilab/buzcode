function [pval, chisq, df, ranks] = FriedmanAnova(M)
  
  ranks = zeros(size(M));
  
  N = size(M, 1);
  k = size(M, 2);
  
  for i = 1:size(ranks, 1)
    ranks(i,:) = rank_order(M(i,:));
  end
  
  chisq = sum(sum(ranks).*sum(ranks)) *12 / (N * k * (k+1)) - (3 * N * (k+1));
%  keyboard
  
  df = k-1;
  
  pval = 1 - chisqp(chisq, df);
  