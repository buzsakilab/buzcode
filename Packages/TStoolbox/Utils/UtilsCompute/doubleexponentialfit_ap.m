function beta=exponentialfit_ap(t,y,beta0)

%  beta=exponentialfit_ap(t,y,beta0)
%  fit to exponential function by calling nlinfit with exponential_ap
%  where f = beta(1)*exp(beta(2)*x) + beta(3);

    beta=nlinfit(t,y,'doubleexponential_ap',beta0);

