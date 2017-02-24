parent_dir = '/home/fpbatta/Data/DIRAC';



datasets = List2Cell([ parent_dir filesep 'dirs_BD1.list' ] );
datasets = List2Cell([ parent_dir filesep 'datasets_mate.list' ] );



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
%A = run(A, 'DIRACRuns', datasets, 'DIRACRuns', 'DoDebug', 1);
%A = run(A, 'ReactRPairs', datasets, 'ReactRateRPairsSPW', 'Debug', 1);
%A = run(A, 'ReactRPairs2', datasets, 'ReactRateRPairs', 'Debug', 1);
%A = run(A, 'ReactRPairs', datasets, 'ReactRateRPairsNoSPW', 'Debug', 1);
%A = run(A, 'DIRACInterneurons', datasets, 'DIRACInterneurons', 'Debug', 1);
%A = run(A, 'DIRACInterneurons2', datasets, 'DIRACInterneurons2', 'Debug', 1);
%A = run(A, 'ReactRatePyrDIRAC', datasets, 'ReactRatePyrDIRAC', 'Debug', 1);
%A = run(A, 'DIRACPrePostSleep', datasets, 'DIRACPrePostSleep', 'Debug',
%1);
%A = run(A, 'DIRACPrePostSleepRates', datasets, 'DIRACPrePostSleepRates', 'Debug', 1);
%A = run(A, 'DIRACDatasets', datasets, 'DIRACDatasets', 'Debug', 1);
%  A = run(A, 'DIRACSimpleMeasuresReact', datasets,
%  'DIRACSimpleMeasuresReact',...
%  	'Debug', 1);
% $$$ A = run(A, 'DIRACSleepACG', datasets, 'DIRACSleepACG', 'Debug', 1);
%A = run(A, 'DIRACPosition', datasets, 'DIRACPosition', 'DoDebug', 1);
%A = run(A, 'DIRACThetaPhase', datasets, 'DIRACThetaPhase', 'DoDebug', 1);
%A = run(A, 'DIRACPhasePlaceCells', datasets, 'DIRACPhasePlaceCells', 'DoDebug', 1);
% A = run(A, 'DIRACTimeCourse', datasets, 'DIRACTimeCourse',...
%   	'DoDebug', 1);
 A = run(A, 'dataToMate2', datasets, 'dataToMate2', 'DoDebug', 1);
 