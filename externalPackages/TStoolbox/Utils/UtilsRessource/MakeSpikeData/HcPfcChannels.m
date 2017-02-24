function A = HcPfcChannels(A)


A = registerResource(A, 'HcChannels', 'numeric', {[], []}, ...
    'hcChannels', ...
    ['Hc eeg channels']);


A = registerResource(A, 'PfcChannels', 'numeric', {[], []}, ...
    'pfcChannels', ...
    ['Hc eeg channels']);


c = current_dataset(A);
c = str2num(c(4:5));

if c==12
	hcChannels = 1;
	pfcChannels = [[2:4] 6]';
else
	hcChannels = [5:6]';
	pfcChannels = [1:4]';
end

A = saveAllResources(A);
