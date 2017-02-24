function [r] = rank_order(x);

u = unique(x);

r = zeros(size(x));

n_used = 0;

for i = 1:length(u)
  ix = find(x == u(i));
  rk = n_used + mean(1:length(ix));
  r(ix) = rk; 
  n_used = n_used + length(ix);
end
