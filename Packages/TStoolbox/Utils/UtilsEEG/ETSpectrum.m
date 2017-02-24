function [ets dt] = ETSpectrum(t1,tspec,spec,bin,nbbins)


nbbins = 2*floor(nbbins/2)+1;
ets = zeros(nbbins,size(spec,2));

for ii=1:size(spec,2)
    [ets(:,ii) dt] = ETAverage(t1,tspec,spec(:,ii),bin,nbbins);
end