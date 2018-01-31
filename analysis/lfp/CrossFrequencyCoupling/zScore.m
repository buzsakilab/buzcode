function zval = zScore(val, mean, sd)
% zval = zScore(val, mean, sd)

zval = (val - mean)./sd;
