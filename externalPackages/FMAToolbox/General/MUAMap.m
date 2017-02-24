function MUAMap(pos,nShanks)

for ii=1:nShanks
  spikesii=GetSpikeTimes([ii -1]);
  [map,stats]=FiringMap(pos,spikesii,'smooth',1,'nBins',[200 100]);
  subplot(7,3,ii)
  PlotColorMap(map.rate,map.time)
  set(gca,'YDir','Reverse')
  clim([0 35])
  xlabel([int2str(ii)]);
end

  