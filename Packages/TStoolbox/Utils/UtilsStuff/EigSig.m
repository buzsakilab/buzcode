function H = EigSig(M)

Ms = M;

M = Ms;
M = zscore(M(1:end,:));
q = size(M,2);
p = size(M,1);
eigV = [];
sigma = [];

for i=1:100

	Mp = [];
	
	for j=1:q
		rp = randperm(p);
		Mp = [Mp M(rp,j)];
	end
	c = corrcoef(Mp);
	c(isnan(c))=0;
%  	eigV = [eigV;eig(c)];
	e = eig(c);
	eigV = [eigV;e/sqrt(mean(e.*e))];

	sigma = [sigma;std(M(~diag(diag(c))))];

end

%  s = 1/sqrt(p);
s = 1;
lambdaMx =p 2*s*sqrt(q); %p=q, lambdaMn = 0;
[C,B] = hist(eigV,100);
rho = real(sqrt((lambdaMx^2-(B-1).^2))./(pi*s));

figure(1),clf
hold on
bar(B,C/length(eigV))
plot(B,rho/sum(rho),'r')

keyboard
