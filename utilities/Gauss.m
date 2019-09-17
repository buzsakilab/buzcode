function [gauss] = Gauss(x,mu,s)
%Gauss(x,mu,s) returns a normalized gaussian with mean mu and standard 
%deviation s

p1 = -.5 * ((x - mu)/s) .^ 2;   %Term in the exponent
p2 = (s * sqrt(2*pi));          %Normalization

gauss = exp(p1) ./ p2;          %Put it all together!

end