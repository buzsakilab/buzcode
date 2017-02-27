function Strengths = SpikeTransfersInPairsPerTimeBin(S,intervals,TotalSeconds,SecondsPerInterval,pairs,funcsynapses)

% SecondsPerInterval = 300;
% TotalSeconds = End(intervals{7});
SampleRate = 10000;%based on TSToolbox

PointsPerInterval = SecondsPerInterval*SampleRate;
TotalPoints = TotalSeconds*SampleRate;


IntervalStarts = [0:PointsPerInterval:TotalPoints-PointsPerInterval];
IntervalEnds = [PointsPerInterval:PointsPerInterval:TotalPoints];
IntervalMeans = mean([IntervalStarts;IntervalEnds],1);
Ints = intervalSet(IntervalStarts,IntervalEnds);

SubS = {};
% subdivide S into time intervals
wh = waitbar(0,'Dividing Spike Trains');
for a = 1:length(length(Ints))
    SubS{a} = Restrict(S,subset(Ints,a));
    Rates(:,a) = Rate(SubS{a});
    Strengths(:,a) = spikeTransfer_InPairSeries(SubS{a},pairs,funcsynapses.BinMs,funcsynapses.CnxnStartTimesVsRefSpk,funcsynapses.CnxnEndTimesVsRefSpk)';
    waitbar(a/length(length(Ints)),wh)
end
close(wh)

% for a = 1:size(pairs,1)
%     avgTransfers(a) = funcsynapses.TransferWeights(pairs(a,1),pairs(a,2));
% end
% avgTransfers = repmat(avgTransfers',[1,size(Strengths,2)]);
% Strengths = Strengths./avgTransfers;
% 

meanStrengths = mean(Strengths,1);
stdStrengths = std(Strengths,[],1);

h = figure;
ax = subplot(2,1,1);
plot(IntervalMeans,meanStrengths);
hold on
plot(IntervalMeans,smooth(meanStrengths',20),'k')
plotIntervalsStrip(ax,intervals)
title('Synaptic Strength Mean')

ax = subplot(2,1,2);
plot(IntervalMeans,stdStrengths);
hold on
plot(IntervalMeans,smooth(stdStrengths',20),'k')
plotIntervalsStrip(ax,intervals)
title('Synaptic Strength StandardDev')
