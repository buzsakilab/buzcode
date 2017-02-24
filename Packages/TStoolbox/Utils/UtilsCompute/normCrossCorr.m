function [H,B] = normCrossCorr(tx,ty,bins,nbBins,epoch)

mx = rate(tx,epoch);
my = rate(ty,epoch);
T = tot_length(epoch,'s');

[H,B] = CrossCorr(Range(tx),Range(ty),bins,nbBins);
H = mx*(H - my);
H = sqrt(bins/1000*T/(mx*my))*H;
