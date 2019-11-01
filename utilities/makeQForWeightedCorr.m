function Q = makeQForWeightedCorr(UniqueNumBins, numPlaces)
%function Q = makeQForWeightedCorr(UniqueNumBins, numPlaces)

a1 = repmat((1:max(UniqueNumBins))', [1, numPlaces]);
b1 = repmat(1:numPlaces, [max(UniqueNumBins), 1]);
Q = {};
for i = 1:length(UniqueNumBins)
    
    Q{UniqueNumBins(i)} = [reshape(a1(1:UniqueNumBins(i), :), [], 1), reshape(b1(1:UniqueNumBins(i), :), [], 1)];
end