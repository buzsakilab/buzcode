function [normX,Bx,By]=normalizepdf(X)
% Forces the pdf of data to have a normal distribution using a data adaptive lookup table. 
% 
% [normX,Bx,By]=normalizepdf(X)
%
% normX=N(X) where N is an data adaptive monotonically increasing function.
% normX will have zero mean and unit variance. 
% 
% Bx,By describes the lookup table
% 
% (c) Aslak Grinsted 2002-2004

[Bx,Js]=sort(X(:)); 
n=size(Bx,1);
d=(diff(Bx)~=0);

I=(1:n)';
I=I([d;true]);
%if (nargout>1)
    Bx=Bx(I);
%end

J=cumsum([1;d]);

n=length(X);
I=[0;I];

By=(.5/n)*(I(1:(end-1))+(I(2:end))); %percentile
By=norminv(By,0,1);

normX=interp1q(Bx,By,X);