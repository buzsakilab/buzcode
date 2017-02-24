function P = CorrelationCompare(r1,r2,N1,N2)

%  P = CorrelationCompare(r1,r2,N1,N2)
%  Statistical comparison of two correlations coefficients.
%  INPUT:
%  	r1,2: correlation coefficients to compare
%  	N1,2: number of samples in population 1,2
%  OUTPUT:
%  	P, the p-vlue that r1 is different from r2;
%  
%  Adrien Peyrache, 2008


Zf1 = 1/2 * log( (1+r1) / (1-r1) );
Zf2 = 1/2 * log( (1+r2) / (1-r2) );

z = (Zf1 - Zf2) / sqrt( 1/(N1-3) + 1/(N2-3) );

P = normpdf(z);
