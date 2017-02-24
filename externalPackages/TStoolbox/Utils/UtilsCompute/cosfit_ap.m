function beta=cosfit_ap(t,y,beta0)

%  beta=cosfit_ap(t,y,beta0)
%  fit to cosine function by calling nlinfit with cos_ap
%try
    beta=nlinfit(t,y,'cos_ap',beta0);
%catch
%    beta = [NaN NaN];
%end
