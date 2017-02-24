function A = PhasePref(A)


A = getResource(A,'SpikeData');
A = getResource(A,'HcChannels');
A = getResource(A,'MazeEpoch');
mazeEpoch =  mazeEpoch{1};

A = registerResource(A, 'HcThetaPhase', 'cell', {[], []}, ...
    'hcThetaPhase', ...
    ['spike phase with hc theta'],'mfile');

A = registerResource(A, 'HcThetaRhythmPhase', 'cell', {[], []}, ...
    'hcThetaRhythmPhase', ...
    ['spike phase with hc theta'],'mfile');


hcThetaPhase = cell(6,1);
hcThetaRhythmPhase = cell(6,1);
[d,ds,d] = fileparts(current_dir(A));

for hcTrace=hcChannels'
	
	eegfname = [current_dir(A) filesep ds 'eeg' num2str(hcTrace) '.mat'];

	try 	
		if exist([eegfname '.gz'])
		display(['unzipping file ' eegfname]);
		eval(['!gunzip ' eegfname '.gz']);
		end
		load(eegfname)
		display(['zipping file ' eegfname]);
		eval(['!gzip ' eegfname]);
		
		eval(['eegHc = Restrict(EEG' num2str(hcTrace) ',mazeEpoch);']);
		eval(['clear EEG' num2str(hcTrace) ';']);
	
		dEegHc = resample(Data(eegHc),1,10);
		rg = Range(eegHc);
		rg = rg(1:10:end);
	
		Fn = 1/(2*median(diff(rg/10000)));
		
		b = fir1(96,[5 10]/Fn);
		dEegHc = filtfilt(b,1,dEegHc);
		
		if length(rg) ~= length(dEegHc)
			keyboard;
		end
		
		eegHc = tsd(rg,dEegHc);
		S = Restrict(S,mazeEpoch);	
	
		[phT phS] = firingPhaseHilbert(eegHc,S);
		hcThetaPhase{hcTrace} = phS;
		hcThetaRhythmPhase{hcTrace} = phT;

	catch

		warning(['problem in ' ds '. TT = ' num2str(hcTrace) '; ' lasterr])
	
	end
end

A = saveAllResources(A);
