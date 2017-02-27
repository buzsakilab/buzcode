function [boxedX,Bx,By]=boxpdf(X)
% Forces the pdf of data to have a boxed distribution using a data adaptive lookup table. 
% 
% [boxedX,Bx,By]=boxpdf(X)
%
% boxedX=N(X) where N is an data adaptive monotonically increasing function.
% boxedX vary between zero and one. 
% 
% Bx,By describes the lookup table
% 
% (c) Aslak Grinsted 2002-2004


%tested with: [prctile(x,boxpdf(x)*100)-x]

%home made fast unique:
[Bx,Js]=sort(X(:)); 
n=size(Bx,1);
d=(diff(Bx)~=0);

I=(1:n)';
I=I([d;logical(1)]);
%if (nargout>1)
    Bx=Bx(I);
%end

J=cumsum([1;d]);

n=length(X);
I=[0;I];

By=(.5/n)*(I(1:(end-1))+(I(2:end))); %percentile


%boxedX=By(J(Js([2:end 1])));
boxedX=interp1q(Bx,By,X);
% 
% %old boxpdf:
% 
% [sx i]=sort(X(:));
% [n, m] = size(sx);
% minx  = min(sx);
% maxx  = max(sx);
% range = maxx-minx;
% 
% 
% % Use the same Y vector if all columns have the same count
% eprob = [0.5./n:1./n:(n - 0.5)./n]';
% 
%  
% %group x's with same value
% I=[1; find(diff(sx)~=0)+1]; %unique indices
% Bx=sx(I); %unique values
% nI=diff([I;length(sx)+1])-1;
% By=zeros(length(Bx),1);
% for i=length(I):-1:1
%     
%     By(i,1)=mean(eprob(I(i)+(0:nI(i))));
% end
% %plot(sx,y);
% 
% 
% 
% boxedX=interp1q(Bx,By,X);
