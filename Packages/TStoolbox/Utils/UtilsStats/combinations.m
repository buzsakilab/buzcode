function b = combinations(n, k)
% b = combinations(n, k)
%
% the binomial coefficient C(n, k)
  
  
if n < k
  error(' n < k ');
end

b = factorial(n) / (factorial(n-k) * factorial(k));