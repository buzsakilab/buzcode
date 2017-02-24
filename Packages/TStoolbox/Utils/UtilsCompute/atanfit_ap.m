function beta=atanfit_ap(t,y,beta0)

%  beta=exponentialfit_ap(t,y,beta0)
%  fit to exponential function by calling nlinfit with exponential_ap
%  where f = beta(1)*exp(beta(2)*x);
try
    beta=nlinfit(t,y,'atan_ap',beta0);
catch
    keyboard
    beta = [NaN NaN];
end
