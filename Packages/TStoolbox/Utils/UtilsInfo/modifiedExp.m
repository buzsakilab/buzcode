function y = modifiedExp(x)

y=exp(x);
y(x>0) = x(x>0)+1;