function [rho,phi] = cart2poln(data)

% [rho,phi] = cart2poln(DATA)
% transforms the p samples of cartesian N-dimensional data (Data is a p-by-n matrix) into spherical coordinates (N-sphere)
% rho is a vector of length p, phi a  matrix of angular values of size p x n-1
% Adrien Peyrache, 2013



squared = data.*data;
N = size(data,2);

rho=sqrt(sum(data,2));
phi = zeros(size(data,1),N-1);
for ii=1:N-1
    phii = data(:,ii)./sqrt(sum(squared(:,ii:end),2));
    phi(:,ii) = acos(phii);
end

% phii = (sqrt(sum(squared(:,end-1:end),2))+data(:,end-1))./data(:,end);
% phi(:,end) = 2*acos(phii);