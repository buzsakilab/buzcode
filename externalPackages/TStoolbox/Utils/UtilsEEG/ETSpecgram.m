function [ets,time,f,S] = ETSpecgram(eeg,evt,nbBins,fq,movingwin);

rg = Range(eeg,'s');
deeg = Data(eeg);
deeg = WhitenSignal(deeg);

%  fq = [55 100];

params.Fs = 1/median(diff(rg));
params.fpass = fq;
params.err = [2, 0.95];
params.trialave = 0;
params.pad = 16;
params.tapers = [3 5];

[S,t,f]=mtspecgramc(deeg,movingwin,params);

t = (t'+rg(1))*10000;

if ~mod(nbBins,2)
nbBins = nbBins+1;
end

ets = zeros(size(S,2),nbBins);

for i=1:size(S,2)

  [eta,dt] = ETAverage(evt,t,S(:,i),movingwin(2)*1000,nbBins);
  ets(i,:) = eta';

end

time = dt;