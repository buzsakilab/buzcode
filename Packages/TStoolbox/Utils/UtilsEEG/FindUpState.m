function A = truc(A)

% Find Upstate 
% Adapted from Muckovski, Cerebral Cortex, 2007
% Adrien Peyrache 2007


A = getResource(A,'PfcTrace');
A = getResource(A,'Sleep1DeltaEpoch');
A = getResource(A,'Sleep1SpindleEpoch');
A = getResource(A,'Sleep2DeltaEpoch');
A = getResource(A,'Sleep2SpindleEpoch');
A = getResource(A,'MidRipS1SWS');
A = getResource(A,'MidRipS2SWS');
A = getResource(A,'SpikeData');
A = getResource(A,'Sleep2Epoch');
A = getResource(A,'MazeEpoch');






%  A = registerResource(A, 'Sleep1UpState', 'cell', {[], []}, ...
%      'sleep1UpState', ['sleep1 sws upstate']);
%  
%  A = registerResource(A, 'Sleep1UpState', 'cell', {[], []}, ...
%      'sleep1UpState', ['sleep1 sws upstate']);

swsEpoch{1} = union(sleep1DeltaEpoch{1},sleep1SpindleEpoch{1});
swsEpoch{2} = union(sleep2DeltaEpoch{1},sleep2SpindleEpoch{1});

[p,dataset,e] = fileparts(current_dir(A));
clear p e;

eegfname = [current_dir(A) filesep dataset 'eeg' num2str(pfcTrace) '.mat'];
if exist([eegfname '.gz'])
	display(['unzipping file ' eegfname]);
	eval(['!gunzip ' eegfname '.gz']);
end
load(eegfname)
display(['zipping file ' eegfname]);
%  eval(['!gzip ' eegfname]);	
	
eval(['eegG = EEG' num2str(pfcTrace)]);
eval(['clear EEG' num2str(pfcTrace)]);

dEeg = resample(Data(eegG),1,5);
rg = Range(eegG);
rg = rg(1:5:end);
eegG = tsd(rg,dEeg);
	
clear dEeg;

%  eeg{1} = Restrict(eegG,swsEpoch{1});
eeg{1} = Restrict(eegG,mazeEpoch{1});
eeg{3} = Restrict(eegG,sleep2Epoch{1});
eeg{2} = Restrict(eegG,subset(swsEpoch{2},1));

clear eegG;

bin = 5; %in ms
win = 100; %#nb of points

Q = MakeQfromS(S,100);

e = subset(swsEpoch{2},1);


for s=1:3
		
		b = fir1(256,[48 52]/200,'stop');
		dEegF50 = filtfilt(b,1,Data(eeg{s}));

		b = fir1(96,[40 100]/200);
		dEegF = filtfilt(b,1,dEegF50);
		
		
		dEeg = dEegF.^2;
		dEeg = convn(dEeg,gausswin(10),'same')/sum(gausswin(10));
		dEegF = sqrt(dEeg);
		dEegF = convn(dEegF,gausswin(100),'same')/sum(gausswin(100));
	
		figure(s),clf
		[C{s},B] = hist(dEegF,100);
		bar(B,C{s})
		m(s)= mean(dEegF);
		sd(s)=std(dEegF);
		title(['mean=' num2str(m(s)) ' ; std=' num2str(sd(s)) ]);

		if 0
	
			eegF = tsd(Range(eeg{s}),dEegF);
			upStG = thresholdIntervals(tsd(Range(eeg{s}),zscore(dEegF)),0.5)
			upStG = mergeCloseIntervals(upStG,5000);
			dwStG = thresholdIntervals(tsd(Range(eeg{s}),zscore(dEegF)),0,'Direction','Below')
			dwStG = mergeCloseIntervals(dwStG,5000);
	
			stUp = Start(upStG);
			stDw = Start(dwStG);
	
			for i = 1:length(stDw)-1
	
				if length(find((stUp>stDw(i) & stUp<stDw(i+1))))==0
					stDw(i+1)=stDw(i);
				end
	
			end
			
			for i = 1:length(stUp)-1
	
				if length(find((stDw>stUp(i) & stDw<stUp(i+1))))==0
					stUp(i+1)=stUp(i);
				end
	
			end
	
			warning off
			upStG = intervalSet(stUp,End(upStG));
			upStG = mergecloseIntervals(upStG,1);
			dwStG = intervalSet(stDw,End(dwStG));
			dwStG = mergecloseIntervals(dwStG,1);
			warning on	
	
			Qr = Restrict(Q,e);
			dQ = mean(zscore(full(Data(Qr)))')';
			dQ = convn(dQ,gausswin(10),'same')/sum(gausswin(10));
					

			figure,clf,hold on
			plot(Range(Qr,'s'),dQ/10+0.1,'r')
			rg = Range(eeg{s},'s');
			plot(rg,Data(eeg{s})/20+0.2,'k')
			plot(Range(eeg{s},'s'),dEegF,'k')
		
			plot(Range(Restrict(midRipS2SWS{1},e),'s'),0.1,'r*')
			line([rg(1) rg(end)],[m(1)-3*sd(1) m(1)-3*sd(1)])
%  			plot(Start(upStG,'s'),-2.7,'m*')
%  			plot(Start(dwStG,'s'),-2.9,'b*')
		
		end

end
	keyboard


%  sleep1UpState = {upState{1}};
%  sleep2UpState = {upState{2}};
%  
%  A = saveAllResources(A);


