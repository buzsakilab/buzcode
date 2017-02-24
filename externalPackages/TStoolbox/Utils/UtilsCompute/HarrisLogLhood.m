function L = HarrisLogLhood(s,Q)

sa = s*w;
n = Q(:,cellIx);


end

function y = g(x)
if x<0
  y = exp(x);
else
  y = x+1;
end