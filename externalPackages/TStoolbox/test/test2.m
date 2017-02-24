parent_dir = '/home/fpbatta/Data/WarpRat';



datasets = List2Cell([ parent_dir filesep 'dirs_rate_react.list' ] );

%datasets = datasets(12:end);

A = Analysis(parent_dir);

%A = run(A, 'SpwClusters', datasets, 'SpwClusters', 'Debug', 1);

A = run(A, 'SpwPETHClusters', datasets, 'SpwPETHClusters', 'Debug', 1);
%A = run(A, 'SpwPETHTemplate', datasets, 'SpwPETHTemplate', 'Debug', 0);

%A = run(A, 'ReactRateSpwClusteredPETH', datasets, 'ReactRateSpwClusteredPETH', 'Debug', 1);
