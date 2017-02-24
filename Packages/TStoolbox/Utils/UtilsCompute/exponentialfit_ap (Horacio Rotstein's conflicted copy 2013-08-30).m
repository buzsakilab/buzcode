function [beta,r,J,Sigma,mse,errorparam,robustw]=exponentialfit_ap(t,y,beta0)
try
    beta=nlinfit(t,y,'exponential_ap',beta0);
catch
    beta = [NaN NaN];
end
