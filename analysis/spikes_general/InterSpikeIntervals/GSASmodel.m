function [allISIdist] = GSASmodel(GSASparms,logtbins,numcells,numAS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%Note: vectr form uses logCV... for fitting
%logtbins should be base e

if ~isstruct(GSASparms)
    GSASparms = convertGSASparms(GSASparms,numcells,numAS);
end

%Here: if logtbins = 'sample'. Put in a very large vector of possible times
%to sample from, save that we're sampling.
sample = false;
if strcmp(logtbins,'sample')
    sample = true;
    logtbins = linspace(-10,7,1000);
end

GSISI = LogGamma(GSASparms.GSlogrates,GSASparms.GSCVs,GSASparms.GSweights,logtbins');
%%
ASISI = zeros([size(GSISI),length(GSASparms.ASlogrates)]);
for aa = 1:length(GSASparms.ASlogrates)
    ASISI(:,:,aa) = LogGamma(GSASparms.ASlogrates(aa),GSASparms.ASCVs(aa),GSASparms.ASweights(:,aa)',logtbins');
end
%%
allISIdist = sum(ASISI,3)+GSISI;

if sample
    numsamps = 30000;
    allISIdist = randsample(logtbins,numsamps,true,allISIdist); 
end

end

