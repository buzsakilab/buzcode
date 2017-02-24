hm = getenv('HOME');


parent_dir =  [hm filesep 'Data/amyphys'];
cd(parent_dir);


datasets = List2Cell([ parent_dir filesep 'datasets_amyphys.list' ] );

%

A = Analysis(parent_dir);


% A = run(A, 'AmyphysCellList_fn', datasets, ...
%       	'AmyphysCellList_fn', 'DoDebug', 1);
% 
%     keyboard
% A = run(A, 'AmygdalaSpikeData_fn', datasets, ...
%  	'AmygdalaSpikeData_fn', 'DoDebug', 1);
%  A = run(A, 'AmyphysDoubles', datasets, 'AmyphysDoubles', 'DoDebug', 1);
%  A = run(A, 'AmyPhysStableInterval', datasets, ...
%      	'AmyPhysStableInterval', 'DoDebug', 1);
% 
%  A = run(A, 'AmyphysDatasets', datasets, ...
%   	'AmyphysDatasets', 'DoDebug', 1);
% 
% A = run(A, 'AmyphysBehavData', datasets, ...
%    	'AmyphysBehavData', 'DoDebug', 1);
%  
%  A = run(A, 'AnatomyData_fn', datasets, ...
%   	'AnatomyData', 'DoDebug', 1);

% A = run(A, 'AmyphysParseTrials', datasets, ...
%   	'AmyphysParseTrials', 'DoDebug', 1);

% A = run(A, 'AmyphysBaselineRateBurst', datasets, ...
%  	'AmyphysBaselineRateBurst', 'DoDebug', 1);
% 
% A = run(A, 'AmyphysVisualResponsiveness', datasets, ...
%  	'AmyphysVisualResponsiveness', 'DoDebug', 1);

%  A  =run(A, 'AmyphysVisualPETH3', datasets, ...
%   	'AmyphysVisualPETH3', 'DoDebug', 1);
% % 
% A  =run(A, 'AmyphysVisualPETH2', datasets, ...
%   	'AmyphysVisualPETH2', 'DoDebug', 1);
% A  =run(A, 'AmyphysAllAnovas', datasets, ...
%      	'AmyphysAllAnovas', 'DoDebug', 1);
% % 
% A = run(A, 'AmyphysVisualInfo', datasets, 'AmyhpysVisualInfo', 'DoDebug',
% 1);

%A = run(A, 'AmyphysStimSet', datasets, 'AmyphysStimSet', 'DoDebug', 1);

%A = run(A, 'AmyphysStimRate', datasets, 'AmyphysStimRate', 'DoDebug', 1);
A = run(A, 'AmyphysDistMatrix', datasets, 'AmyphysDistMatrix', 'DoDebug', 1);
%A = run(A, 'AmyphysCatSpec', datasets, 'AmyphysCatSpec', 'DoDebug', 1);

%A = run(A, 'AmyphysExampleHist', datasets, 'AmyphysExamples', 'DoDebug',1);
% A = run(A, 'AmyphysMonkeyStimRate', datasets, 'AmyphysMonkeyStimRate', 'DoDebug', 1);
% 
% A = run(A, 'AmyphysSelectAnovas', datasets, 'AmyphysSelectAnovas', 'DoDebug', 1);