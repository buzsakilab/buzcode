function A = truc(A)

% Find Upstate 
% Adapted from Muckovski, Cerebral Cortex, 2007
% Adrien Peyrache 2007


A = getResource(A,'PfcTrace');
A = getResource(A,'Sleep1DeltaEpoch');
A = getResource(A,'Sleep1SpindleEpoch');
A = getResource(A,'Sleep2DeltaEpoch');
A = getResource(A,'Sleep2SpindleEpoch');

load([current_dir(A) filesep 'SleepLarge.mat;']);

s2 = sleepLargeSpecgram{2};
d = Data(s2);
rg = Range(s2);
hf = find(sleepLargeSpecgramFreq>65);

hp = mean(log10(d(:,hf)'+eps));

