function f=autocorr_index(p,x)

% f=autocorr_index(params,x)
% f=(p(1)*(sin(2*pi*p(4)*x)+1)+p(2))*exp(-abs(x)/p(5))+p(3)*exp(-abs(x).^2/p(5)^2);
% theta index: p(1)/p(2)
% c.f. Royer et al., 2010

f=(p(1)*0.5*(cos(2*pi*p(4).*x)+1)+p(2)).*exp(-abs(x)/p(5))+p(3).*exp(-abs(x).^2/p(6)^2);
%f=(p(1)*(cos(2*pi*p(3).*x)+1)+p(2)).*exp(-abs(x)/p(4));

