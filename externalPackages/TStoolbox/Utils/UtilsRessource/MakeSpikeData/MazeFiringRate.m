function A = truc(A)

A = getResource(A,'MazeEpoch');
A = getResource(A,'Sleep1Epoch');
A = getResource(A,'Sleep2Epoch');
A = getResource(A,'SpikeData');

A = registerResource(A, 'MeanFr', 'numeric', {[], []}, ...
    'meanFr', ['mean firing rate']);

meanFr = [];
try
	ep = union(sleep1Epoch{1},mazeEpoch{1},sleep2Epoch{1});
catch
	keyboard
end


for i=1:length(S)

	meanFr = [meanFr;mean(Data(intervalRate(S{i},ep)))];

end

A = saveAllResources(A);