function [pval, chisq, df] = CircularANOVA(th)
% [pval, chisq, df] = CircularANOVA(sr) Circular Non-parametric ANOVA
%
% INPUTS:
% th: a cell array, each for each group, each being an array of angular
% values
% OUTPUTS 
% pval: the p-value of the ANOVA
% chisq: the chi-squared value
% df: the number of degrees of freedom 


for i = 1:length(th)
    n(i) = length(th{i});
    [mu(i), rm, delta(i), p] = CircularMean(th{i});
end

if max(delta) / min(delta) < 4 % Fisher p. 116
    display('Using method 1');
    Cp = sum(n .* cos(mu));
    Sp = sum(n.* sin(mu));
    Rp = sqrt(Cp^2 + Sp^2);
    N = sum(n);
    delta0 = sum(n .* delta) / N;
    
    chisq = 2 * (N-Rp)/delta0;
    
    
else
    display('Using method 2');
    sigma2 = delta ./ n;
    
    Cm = sum(cos(mu) ./ sigma2);
    Sm = sum(sin(mu) ./ sigma2);
    Rm = sqrt(Cm^2 + Sm^2);
    
    chisq = 2 * (sum(1./sigma2) - Rm);
    

end


df = length(th)-1;

pval = 1 - chi2cdf(chisq, df);


