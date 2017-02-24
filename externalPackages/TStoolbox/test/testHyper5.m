parent_dir = '/home/fpbatta/Data/Hyper5';



datasets = List2Cell([ parent_dir filesep 'dirs_Hyper5.list' ] );

%datasets = datasets(12:end);

A = Analysis(parent_dir);

%A = run(A, 'SpwClusters', datasets, 'SpwClusters', 'Debug', 1);

%A = run(A, 'DIRACCellList_fn', datasets, 'DIRACCellList', 'Debug', 1);
%A = run(A, 'HippoSpikeData_BD1_fn', datasets, 'SpikeData', 'Debug', 1);
%A = run(A, 'SpwPETHTemplate', datasets, 'SpwPETHTemplate', 'Debug', 0);
%A = run(A, 'DIRACEpochs10min', datasets, 'DIRACEpochs10min', 'Debug', 1);
%A = run(A, 'ReactRateSpwClusteredPETH', datasets, 'ReactRateSpwClusteredPETH', 'Debug', 1);

%A = run(A, 'DIRACEpochs10min', datasets, 'DIRACEpochs10min', 'Debug', 1);
%A = run(A, 'DIRACSpw', datasets, 'DIRACSpw', 'Debug', 1);
%A = run(A, 'ReactRateDIRAC', datasets, 'ReactRateDIRAC', 'Debug', 1);
%A = run(A, 'DIRACRuns', datasets, 'DIRACRuns', 'Debug', 1);
%A = run(A, 'ReactRPairs', datasets, 'ReactRateRPairs', 'Debug', 1);
%A = run(A, 'DIRACInterneurons', datasets, 'DIRACInterneurons', 'Debug', 1);
%A = run(A, 'DIRACInterneurons2', datasets, 'DIRACInterneurons2', 'Debug', 1);
%A = run(A, 'ReactRatePyrDIRAC', datasets, 'ReactRatePyrDIRAC', 'Debug', 1);
%A = run(A, 'Hyper5SpikeData_fn', datasets, 'Hyper5SpikeData', 'Debug', 1);

%A = run(A, 'Hyper5Epochs10min', datasets, 'Hyper5Epochs10min', 'Debug', 1);
%A = run(A, 'Hyper5Spw', datasets, 'Hyper5Spw', 'Debug', 1);
%A = run(A, 'ReactRateHyper5', datasets, 'ReactRateHyper5', 'Debug', 1);
%A = run(A, 'Hyper5PrePostSleep', datasets, 'Hyper5PrePostSleep',
%'Debug', 1);
%A = run(A, 'ReactRPairsHyper5', datasets, 'ReactRPairsPyrHyper5', 'Debug', 1);
%A = run(A, 'Hyper5Datasets', datasets, 'Hyper5Datasets', 'Debug', 1);
% $$$ A = run(A, 'Hyper5SimpleMeasuresReact', datasets, ...
% $$$ 	'Hyper5SimpleMeasuresReact', 'Debug', 1);

% A = run(A, 'Hyper5ArmTimes', datasets, 'Hyper5ArmTimes', ...
% 	'Debug', 1);
A = run(A, 'Hyper5TimeCourse', datasets, 'Hyper5TimeCourse', ...
	'DoDebug', 1);
