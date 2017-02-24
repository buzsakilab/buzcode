function beta=exponentialfit_ap(t,y,beta0)

%  beta=exponentialfit_ap(t,y,beta0)
%  fit to exponential function by calling nlinfit with exponential_ap
%  where f = beta(1)*exp(beta(2)*x) + beta(3);
try
    beta=nlinfit(t,y,'exponential_ap',beta0);
catch
    beta = [NaN NaN NaN];
end
