function A = truc(A)

A = getResource(A,'HcTraceRipples');
A = getResource(A,'PfcTrace');

A = getResource(A,'Sleep1Epoch');
A = getResource(A,'Sleep2Epoch');

A = registerResource(A,'Sleep1RipplesRMS', 'tsdArray', {[], []}, 'sleep1RipplesRMS', 'ripples RMS');
A = registerResource(A,'Sleep2RipplesRMS', 'tsdArray', {[], []}, 'sleep2RipplesRMS', 'ripples RMS');

[p,dataset,e] = fileparts(current_dir(A));
clear p e;

eegfname = [current_dir(A) filesep dataset 'eeg' num2str(hcTraceRipples) '.mat'];
if exist([eegfname '.gz'])
	display(['unzipping file ' eegfname]);
	eval(['!gunzip ' eegfname '.gz']);
end
load(eegfname)
display(['zipping file ' eegfname]);
eval(['!gzip ' eegfname]);	
	
eval(['eegG = EEG' num2str(hcTraceRipples)]);
eval(['clear EEG' num2str(hcTraceRipples)]);

for s=1:2
	eval(['eeg = Restrict(eegG,sleep' num2str(s) 'Epoch{1});']);
	dEeg = Data(eeg);
	rg = Range(eeg);
	
	b = fir1(96,[0.1 0.3]);
	dEeg = filtfilt(b,1,dEeg);
	dEeg = dEeg.^2;
	dEeg = convn(dEeg,gausswin(20),'same')/sum(gausswin(20));
	eval(['sleep' num2str(s) 'RipplesRMS = {tsd(rg,sqrt(dEeg))};']);
end

A = saveAllResources(A);

