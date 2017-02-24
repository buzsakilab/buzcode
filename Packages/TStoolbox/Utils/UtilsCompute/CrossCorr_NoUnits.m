function [C,B] = CrossCorr(t1,t2,binsize,nbins)


%[C,B] = CrossCorr_NoUnits(t1,t2,binsize,nbins)
%
% a timeless point process crosscorr function without any assumption on time units,
% INPUTS:
%     t1,2: vectors of event timing
%     binsize: width of time bin in the same units as t1,2
%     nins: total number of bins
% OUTPUT:
%     C: frequency of occurence of events in t2 relative to events in t1.
%     The frequency is expressed in the inverse of the time units of t1,2
%     and binsize.
%     B: vector of time-lags.
% 
% Adrien Peyrache


nt1  = length(t1);
nt2 = length(t2);


% we want nbins to be odd */
if floor(nbins / 2)*2 == nbins
    nbins = nbins+1;
end
  
m = - binsize * (nbins / 2);
B = zeros(nbins,1);      
for j = 1:nbins
	B(j) = m + j * binsize;
end
    
% cross correlations */
  
w = (nbins / 2) * binsize;
C = zeros(nbins,1);
i2 = 2;

for i1 = 1:nt1
    lbound = t1(i1) - w;
    while t2(i2) < lbound && i2 < nt2-1
        i2=i2+1;
    end
    while t2(i2-1) > lbound && i2 > 2
	i2=i2-1;
    end
    
    rbound = lbound;
    l = i2;
    for j = 1:nbins
        k = 0;
        rbound = rbound+binsize;
        while t2(l) < rbound && l < nt2-1  
            l = l+1;
            k = k+1;
        end
      
	  C(j) = C(j)+k;
    end
  
end


  for j = 1:nbins
    C(j) = C(j)/(nt1 * binsize);
 
  end
  
