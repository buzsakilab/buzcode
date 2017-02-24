function A = SpindlePhase(A)

%Parameters
b = fir1(96,[0.05 0.1]);
%  cellId = 44;

A = getResource(A,'SpikeData');
A = getResource(A,'PfcTrace');
A = getResource(A,'Sleep1FineSpindleEpoch');
A = getResource(A,'Sleep2FineSpindleEpoch');
A = getResource(A,'Sleep1Epoch');
A = getResource(A,'Sleep2Epoch');
A = getResource(A,'Sleep1DwStE');
A = getResource(A,'Sleep2DwStE');
A = getResource(A,'Sleep1DwStS');
A = getResource(A,'Sleep2DwStS');


A = registerResource(A, 'SpindlesPhase', 'cell', {[], []}, ...
    'spindlesPhase', ...
    ['spike phase with hc theta'],'mfile');


[p,ds,e] = fileparts(current_dir(A));

eegfname = [current_dir(A) filesep ds 'eeg' num2str(pfcTrace) '.mat'];

if exist([eegfname '.gz'])
display(['unzipping file ' eegfname]);
eval(['!gunzip ' eegfname '.gz']);
end
load(eegfname)
%  display(['zipping file ' eegfname]);
%  eval(['!gzip ' eegfname]);

eval(['eegG = EEG' num2str(pfcTrace)]);
eval(['clear EEG' num2str(pfcTrace) ';']);
	
dEeg = Data(eegG);
dEeg = resample(dEeg,1,10);
rg = Range(eegG);
rg = rg(1:10:end);
eegG = tsd(rg,dEeg);
clear dEeg;


eeg{1} = Restrict(eegG,sleep1Epoch{1});
eeg{2} = Restrict(eegG,sleep2Epoch{1});

stDw{1} = Range(sleep1DwStS{1});
stDw{2} = Range(sleep2DwStS{1});

enDw{1} = Range(sleep1DwStE{1});
enDw{2} = Range(sleep2DwStE{1});

spinEp{1} = sleep1FineSpindleEpoch{1};
spinEp{2} = sleep2FineSpindleEpoch{1};

for s=1:2
	

	b = fir1(96,[0.1 0.2]);
	dEegF = filtfilt(b,1,Data(eeg{s}));
	eegF = tsd(Range(eeg{s}),dEegF);
	
	dwEp = intervalSet(stDw{s}-1000,enDw{s}+1000);
	
	spEp = spinEp{s}-dwEp;

	phSpin = ThetaPhase(S, eegF, Start(spEp), End(spEp));
	spindlesPhase{s} = phSpin;

end

spindlesPhase = {spindlesPhase};

A = saveAllResources(A);

