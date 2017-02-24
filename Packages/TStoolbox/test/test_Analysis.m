parent_dir = '/home/fpbatta/Data/WarpRat';



datasets = List2Cell([ parent_dir filesep 'dirs_rate_react.list' ] );

%A = Analysis(parent_dir, 'CortCellDepth', datasets, 'CortCellDepth');
%A = Analysis(parent_dir, 'ReactRate', datasets, 'ReactRate');
% A = Analysis(parent_dir, 'WarpEpochs10min', datasets, 'WarpEpochs10min');
%A = Analysis(parent_dir, 'SpikeData', datasets, 'CortSpikeTrains');
% A = Analysis(parent_dir, 'CortCellList_fn', datasets, 'CortCellList');

%A = Analysis(parent_dir, 'SQSpecResources', datasets, 'SQResources');

A = Analysis(parent_dir);
%A = run(A, 'ReactRateByDelta', datasets, 'ReactRateByDelta');
%A = run(A, 'globalBinnedFRateSleep2', datasets, ...
%	'globalBinnedFRateSleep2');

%A = run(A, 'SleepBurstiness', datasets, 'SleepBurstiness');
%A = run(A, 'showPopSpecgram', datasets, 'SleepBurstiness');
%A = run(A, 'ReactRateSpindles', datasets, 'ReactRateSpindles');

% $$$  A = run(A, 'globalBinnedFRateSleep2', datasets, ...
% $$$  	'globalBinnedFRateSleep2');
%A = run(A, 'ReactRate', datasets, 'ReactRate');
% $$$ A = run(A, 'ReactRateSpwPETH', datasets, 'ReactRateSpwPETH', ...
% $$$ 'Debug', 0);
A = run(A, 'SpwPETH', datasets, 'SpwPETH', 'Debug', 1);

%A = run(A, 'SpwGlobalPETH', datasets, 'SpwGlobalPETH', 'Debug', 1);
