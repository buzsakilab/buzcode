function [pval, chisq, df] = contingency_chisq(M)


[nrow, ncol] = size(M);


sum_c = sum(M, 1);
sum_c = sum_c / sum(sum_c);

sum_r = sum(M, 2);
sum_r = sum_r / sum(sum_r);

exp_freq = sum_r * sum_c * sum(sum(M));

%Mnorm = M / sum(sum(M));

chisq = sum(sum( ((M-exp_freq).^2) ./ exp_freq));

df = (nrow-1)*(ncol-1);

pval = 1 - chisqp(chisq,df);
