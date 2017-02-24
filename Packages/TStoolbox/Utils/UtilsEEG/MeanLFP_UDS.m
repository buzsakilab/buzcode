function A = truc(A)



A = registerResource(A,'PfcLFP_SlowTA','cell',{[], []},'pfcLFP_SlowTA',['pfc lfp trigger averages on UDS']);

A = registerResource(A,'PfcLFP_SlowDeltaTA','cell',{[], []},'pfcLFP_SlowDeltaTA',['pfc lfp trigger averages on UDS']);

A = getResource(A,'PfcChannels');


A = getResource(A, 'Sleep1Epoch');
sleep1Epoch = sleep1Epoch{1};
A = getResource(A, 'Sleep2Epoch');
sleep2Epoch = sleep2Epoch{1};

A = getResource(A, 'Sleep1SlowPeaks');
s1S = sleep1SlowPeaks;
A = getResource(A, 'Sleep2SlowPeaks');
s2S = sleep2SlowPeaks;
A = getResource(A, 'Sleep1SlowDeltaPeaks');
s1SD = sleep1SlowDeltaPeaks;
A = getResource(A, 'Sleep2SlowDeltaPeaks');
s2SD = sleep2SlowDeltaPeaks;

[p,dataset,e] = fileparts(current_dir(A));
clear p e;

minZ = [-2.5:0.25:-0.75 0.75:0.25:2.5];


for ch=1:6

	if sum(ismember(pfcChannels,ch))

		rgS = [Range(s1S{ch});Range(s2S{ch})];
		dS = [Data(s1S{ch});Data(s2S{ch})];
		
		rgSD = [Range(s1SD{ch});Range(s2SD{ch})];
		dSD = [Data(s1SD{ch});Data(s2SD{ch})];

		eegfname = [current_dir(A) filesep dataset 'eeg' num2str(ch) '.mat'];
		if exist([eegfname '.gz'])
			display(['unzipping file ' eegfname]);
			eval(['!gunzip ' eegfname '.gz']);
		end
		load(eegfname)
		display(['zipping file ' eegfname]);
		eval(['!gzip ' eegfname ' &']);	
			
		eval(['eegG = EEG' num2str(ch)]);
		eval(['clear EEG' num2str(ch)]);

		eegG = cat(Restrict(eegG,sleep1Epoch),Restrict(eegG,sleep2Epoch));

		deeg = resample(Data(eegG),1,50);
		rg = Range(eegG);
		rg = rg(1:50:end);

		Fn = 1/(2*median(diff(rg/10000)));
		b = fir1(96,[0.1 4]/Fn);
		deeg = filtfilt(b,1,deeg);
		eegG = tsd(rg,deeg/std(deeg));

		rint = regular_interval(-20000,20000,10000/(Fn));


		for i=1:length(minZ)
			
			display(minZ(i))
			eegTA_S{i,1} = minZ(i);
			eegTA_SD{i,1} = minZ(i);

			if minZ(i)<0
				rgs = rgS(dS<minZ(i));
				rgsd = rgSD(dSD<minZ(i));
			else 
				rgs = rgS(dS>minZ(i));
				rgsd = rgSD(dSD>minZ(i));
			end

			if length(rgs)>0
	
				is = intervalSet(rgs-20000,rgs+20000);
				eegS = oneSeries(intervalSplit(eegG,is,'OffsetStart',-20000));
				eTA = intervalMean(eegS,rint);
				eegTA_S{i,2} = eTA;

			else
				eegTA_S{i,2} = tsd();
			end		

			if length(rgsd)>0

				is = intervalSet(rgsd-20000,rgsd+20000);
				eegS = oneSeries(intervalSplit(eegG,is,'OffsetStart',-20000));
				eTA = intervalMean(eegS,rint);
				eegTA_SD{i,2} = eTA;

			else
				eegTA_SD{i,2} = tsd();
			end

		end

		pfcLFP_SlowTA{ch} = eegTA_S;
		pfcLFP_SlowDeltaTA{ch} = eegTA_SD;

	else

		pfcLFP_SlowTA{ch} = {};
		pfcLFP_SlowDeltaTA{ch} = {};

	end


end

A = saveAllResources(A);