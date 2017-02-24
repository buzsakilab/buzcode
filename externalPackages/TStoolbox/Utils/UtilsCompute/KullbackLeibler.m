function d = KullbackLeibler(c2,c1) 

%  D(P||Q) = KL(Cp,Cq) is the Kullback-Leibler divergence between distribution of covariance Cp and Cq
%  mean are assumed to be 0
%  Adrien Peyrache 2007


N = size(c1,1);
warning off
d = 0.5*(log(det(c2)/det(c1))+trace(inv(c2)*c1) - N);
warning on