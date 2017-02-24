function A = PhasePref(A)

%Parameters
b = fir1(96,[0.05 0.1]);
%  cellId = 44;

A = getResource(A,'CellNames');
A = getResource(A,'SpikeData');
A = getResource(A,'Tetrode');
A = getResource(A,'HcTrace');
%  A = getResource(A,'PfcTrace');
A = getResource(A,'MazeEpoch');
mazeEpoch =  mazeEpoch{1};
%  
%  A = registerResource(A, 'PfcThetaPhase', 'tsdArray', {[], []}, ...
%      'pfcThetaPhase', ...
%      ['spike phase with pfc theta'],'mfile');

A = registerResource(A, 'HcThetaPhase', 'cell', {[], []}, ...
    'hcThetaPhase', ...
    ['spike phase with hc theta'],'mfile');


resdir = [parent_dir(A), filesep 'PhasePlot'];
[p,ds,e] = fileparts(current_dir(A));

hcTrace = {};

for hcTrace=5:6
	
	eegfname = [current_dir(A) filesep ds 'eeg' num2str(hcTrace) '.mat'];
	
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
	dEegHc = filtfilt(b,1,dEegHc);
	
	rg = Range(eegHc);
	rg = rg(1:10:end);
	
	if length(rg) ~= length(dEegHc)
		keyboard;
	end
	
	eegHc = tsd(rg,dEegHc);
	
	thStHc = Start(mazeEpoch);
	thEnHc = End(mazeEpoch);
	
	
	phHc = ThetaPhase(S, eegHc, thStHc, thEnHc);
	hcThetaPhase{hcTrace-4} = phHc;
end

if 0
	
	for i = 1:length(S)
		
		if length(Data(phHc{i}))>0
	
			fh1 = figure(1),clf
			rose(2*pi*Data(phPfc{i}));
			title('Phase Preference with pfc Theta')
			fh2 = figure(2),clf
			rose(2*pi*Data(phHc{i}));
			title('Phase Preference with hc Theta')
			
			saveas(fh1,[resdir filesep ds '_' cellnames{i} 'PhasePrefPfc'],'png');
			saveas(fh2,[resdir filesep ds '_' cellnames{i} 'PhasePrefHc'],'png');
	
		end
	
	end

end


A = saveAllResources(A);
