function [p,chi2] = chi2Table(data)


%  Chi2 test for consitency in a frequency table in 2 conditions (2 class)
%  data is a table in the form sample x class

ni = sum(data')';
n = sum(ni);

xi = data(:,1);
x = sum(xi);

nz = find(ni~=0);
ni = ni(nz);
xi = xi(nz);

K = size(data,1);

chi2 = n^2/(x*(n-x))*(sum((xi.^2)./ni) - x^2/n);
p = chi2pdf(chi2,K-1);