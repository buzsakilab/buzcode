function [eta,etaSem,B] = ETAverage(t1,t2,data,binsize,nbins)

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
    
binsize = binsize*10;  

% cross correlations */
  
w = (nbins / 2) * binsize;
etaSamples = cell(nbins,1);
i2 = 2;

for i1 = 1:nt1
    lbound = t1(i1) - w;
    while t2(i2) < lbound && i2 < nt2
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
            k = k+1;
            l = l+1;
        end
      
	   etaSamples{j} = [etaSamples{j};data(l-k:l-1)];
    end
  
end

eta = zeros(nbins,1);
etaSem = zeros(nbins,1);

for j = 1:nbins
    e = etaSamples{j};
    eta(j) = mean(e);
    etaSem(j) = sem(e);
end
  
