function nSpikes = MClust_CountSpikes(fn)

% MClust_LoadNeuralData
%
%
%  fn = file name string
%
% Only forwards the correct number of arguments
% ADR 2002
% Version M3.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

global MClust_Directory
global MClust_NeuralLoadingFunction

% added ncst 20 May 02
if isempty(MClust_Directory)
    MClustDir = which('MClust.m');
    [MClust_Directory n e] = fileparts(MClustDir);
end

LoadingEngineDirectory = fullfile(MClust_Directory, 'LoadingEngines');
LoadingEngine = MClust_NeuralLoadingFunction;

switch exist(LoadingEngine) %#ok<EXIST>
case {0,1,4,5,7}
    pushdir(LoadingEngineDirectory);
    LoadingEngine = uigetfile( ...
       {'Load*.dll';'Load*.m'}, ...
        'Select a loading engine.');
    if LoadingEngine==0
        skip = 1;
    else
        skip = 0;
    end
    popdir;    
otherwise 
    skip = 0;
end

if skip
    MClust_NeuralLoadingFunction = [];    
    nSpikes = nan;
else
    [LoadingEngineDirectory, LoadingEngine, ext] = fileparts(LoadingEngine);
    MClust_NeuralLoadingFunction = LoadingEngine; 
    nSpikes = feval(LoadingEngine, fn, [], 5);
end
    