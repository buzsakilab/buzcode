parent_dir = '/home/fpbatta/Data/CRAM';



datasets = List2Cell([ parent_dir filesep 'dirs_CRAM.list' ] );

%datasets = datasets(12:end);

A = Analysis(parent_dir);

%A = run(A, 'SpwClusters', datasets, 'SpwClusters', 'Debug', 1);

%A = run(A, 'DIRACCellList_fn', datasets, 'DIRACCellList', 'Debug', 1);
%A = run(A, 'HippoSpikeData_fn', datasets, 'HippoSpikeData', 'Debug', 1);
%A = run(A, 'SpwPETHTemplate', datasets, 'SpwPETHTemplate', 'Debug', 0);
%A = run(A, 'DIRACEpochs10min', datasets, 'DIRACEpochs10min', 'Debug', 1);
%A = run(A, 'ReactRateSpwClusteredPETH', datasets, 'ReactRateSpwClusteredPETH', 'Debug', 1);

%A = run(A, 'DIRACEpochs10min', datasets, 'DIRACEpochs10min', 'Debug', 1);
%A = run(A, 'DIRACSpw', datasets, 'DIRACSpw', 'Debug', 1);
%A = run(A, 'ReactRateDIRAC', datasets, 'ReactRateDIRAC', 'Debug', 1);
%A = run(A, 'DIRACRuns', datasets, 'DIRACRuns', 'Debug', 1);
%A = run(A, 'ReactRPairsCRAM', datasets, 'ReactRateRPairsCRAM', 'Debug', 1);



%A = run(A, 'CRAMPrePostSleep', datasets, 'CRAMPrePostSleep', 'Debug', 1);
%A = run(A, 'CRAMPrePostSleepRates', datasets, 'CRAMPrePostSleepRates', 'Debug', 1);
%A = run(A, 'ReactRPairsCRAM', datasets, 'ReactRateRPairsPyrCRAM', 'Debug', 1);

%A = run(A, 'DIRACInterneurons', datasets, 'DIRACInterneurons', 'Debug',
%1);

%A = run(A, 'CRAMCellList_fn', datasets, 'CRAMCellList', 'Debug', 1);

%A = run(A, 'HippoSpikeData_fn', datasets, 'HippoSpikeData', 'Debug', 1);
%A = run(A, 'DIRACEpochs10min', datasets, 'DIRACEpochs10min', 'Debug',1);
% A = run(A, 'DIRACRuns', datasets, 'DIRACRuns', 'Debug', 1);
%A = run(A, 'ReactRateCRAM', datasets, 'ReactRateCRAM', 'Debug', 1);
%A = run(A, 'ReactRatePyrCRAM', datasets, 'ReactRatePyrCRAM', 'Debug',
%1);

%A = run(A, 'CRAMDatasets', datasets, 'CRAMDatasets', 'Debug', 1);
% A = run(A, 'CRAMSimpleMeasuresReact', datasets, 'CRAMSimpleMeasuresReact', ...
% 	'Debug', 1);

A = run(A, 'CRAMTimeCourse', datasets, 'CRAMTimeCourse', ...
	'DoDebug', 1);

