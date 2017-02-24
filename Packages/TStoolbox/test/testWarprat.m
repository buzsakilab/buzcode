parent_dir = '/home/fpbatta/Data/WarpRat';

datasets = List2Cell([ parent_dir filesep 'place_dir.list' ] );

A = Analysis(parent_dir);

A = run(A, 'WarpPositionAnalysis', datasets, 'WarpPositionAnalysis', 'DoDebug', 1);

