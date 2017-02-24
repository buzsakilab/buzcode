function A = AutoCorrMaze(A)


%Parameters
b = fir1(96,[0.05 0.1]);


A = getResource(A,'MazeEpoch');
mazeEpoch = mazeEpoch{1};

A = getResource(A,'SpikeData');
nbCells = length(S);

A = getResource(A,'CellNames');
A = getResource(A,'HcTrace');

resdir = [parent_dir(A), filesep 'XCorr' filesep 'Maze' filesep 'AutoTheta'];
[p,ds,e] = fileparts(current_dir(A));

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

thresholdHc = median(abs(dEegHc));

[thStHc thEnHc] = findTheta(eegHc,3*thresholdHc,thresholdHc);


for i=1:nbCells

	t = Range(Restrict(S{i},intervalSet(thStHc,thEnHc)));
	[C,B] = CrossCorr(t,t,1,600);
	C(301) = 0;
	
	fh = figure(1),clf
	bar(B,C);
	
	saveas(fh,[resdir filesep ds '_' cellnames{i} 'AutoCorrMazeTheta300ms_1msBin'],'png')

end
