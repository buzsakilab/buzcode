function P = SlopeTest(X1,Y1,X2,Y2);

X1 = X1(:);
X2 = X2(:);

Y1 = Y1(:);
Y2 = Y2(:);

p1 = polyfit(X1,Y1,1);
p2 = polyfit(X2,Y2,1);

y1 = p1(1)*X1+p1(2);
y2 = p2(1)*X2+p2(2);

n1 = length(Y1);
n2 = length(Y2);

Seb1 = sqrt(sum((Y1-y1).^2)/(n1-2))/sqrt(sum((X1-mean(X1)).^2));
Seb2 = sqrt(sum((Y2-y2).^2)/(n2-2))/sqrt(sum((X2-mean(X2)).^2));

Seb12 = sqrt(Seb1^2+Seb2^2);
t = abs(p1(1)-p2(1))/Seb12;

keyboard

if t~=0
%  	P = chi2pdf(t,n1+n2-4);
	P = tpdf(t,n1+n2-4);
else
	P = 1;
end
%  keyboard





