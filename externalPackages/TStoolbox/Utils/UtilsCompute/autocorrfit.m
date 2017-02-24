function beta=autocorrfit(t,y,beta0,betaL,betaU)

%  beta=autocorrfit(t,y,beta0)
% f=(p(1)*(sin(2*pi*p(4)*x)+1)+p(2))*exp(-abs(x)/p(5))+p(3)*exp(-abs(x).^2/p(5)^2);
% theta index: p(1)/p(2)
% c.f. Royer et al., 2010

try
    warning off
    %beta=lsqnonlin(t,y,'autocorr_index',beta0);
    beta = lsqnonlin(@(x) (autocorr_index(x,t)-y),beta0,betaL,betaU);
    
    warning on
catch
    beta = NaN(1,6);
end
