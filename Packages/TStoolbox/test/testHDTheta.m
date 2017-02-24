hm = getenv('HOME');


parent_dir =  [hm filesep 'Data/Angelo'];
cd(parent_dir);


%datasets = List2Cell([ parent_dir filesep 'all_datasets.list' ] );
datasets = List2Cell([ parent_dir filesep 'datasets_eeg.list' ] );
%datasets = datasets(5);

A = Analysis(parent_dir);


%  A = run(A, 'HDThetaCellList_fn', datasets, ...
%    	'HDThetaCellList_fn', 'DoDebug', 1);
% 
%  A = run(A, 'HDThetaSpikeData_fn', datasets, ...
%     	'HDThetaSpikeData', 'DoDebug', 1);
% 
% A = run(A, 'HDThetaPosition', datasets, ...
%    	'HDThetaPosition', 'DoDebug', 1);

% A = run(A, 'HDThetaEEG', datasets, ...
%    	'HDThetaEEG', 'DoDebug', 1);

% A = run(A, 'HDThetaPositionMeasures', datasets, ...
%    	 'HDThetaPositionMeasures', 'DoDebug', 1);
% 


% A = run(A, 'HDThetaThetaEEG', datasets, ...
%    	 'HDThetaThetaEEG', 'DoDebug', 1);
% 

A = run(A, 'HDThetaCellsHDtheta', datasets, ...
    	 'HDThetaCellsHDtheta', 'DoDebug', 1);

% A = run(A, 'HDFlags', datasets, ...
%     	 'HDFlags', 'DoDebug', 1);


% A = run(A, 'HDThetaRatePrecession', datasets, ...
%     	 'HDThetaRatePrecession', 'DoDebug', 1);
%
% A = run(A, 'HDRatePrecFigs', datasets, ...
%     	 'HDRatePrecFigs', 'DoDebug', 1);
%