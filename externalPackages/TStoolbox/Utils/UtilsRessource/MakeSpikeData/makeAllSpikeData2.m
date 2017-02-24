dataDir = '/media/sdc6/Data';

cd(dataDir);

A = Analysis(pwd);

datasets = List2Cell('datasets_noYMPpb.list');
error = []
%  run(A,'makeSpikeData',datasets,'makeSpikeData','DoDebug',0,'Overwrite',0);

for i=2:length(datasets)

	dset = datasets{i};
	cd([dataDir filesep dset])
	
	if ~ exist([dset(7:end) 'eeg1.mat'])
	
		eval(['! son2Mat
		
	else
		error = [error;dset]

	end
end

