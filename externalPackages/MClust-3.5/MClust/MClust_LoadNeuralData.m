function [T,WV] = MClust_LoadNeuralData(fn, records_to_get, key)

% MClust_LoadNeuralData
%
%
%  fn = file name string
%  records_to_get = an array that is either a range of values
%  key = 1: timestamp list.(a vector of timestamps to load 
%               (uses binsearch to find them in file))
%        2: record number list (a vector of records to load)
%        3: range of timestamps (a vector with 2 elements: 
%                a start and an end timestamp)
%        4: range of records (a vector with 2 elements: 
%                a start and an end record number)
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

LoadingEngine = MClust_NeuralLoadingFunction;

switch exist(LoadingEngine) %#ok<EXIST>
case {0,1,4,5,7}
    LoadingEngineDirectory = fullfile(MClust_Directory, 'LoadingEngines');
    pushdir(LoadingEngineDirectory);
    LoadingEngine = uigetfile( ...
       {'Load*.mex';'Load*.dll';'Load*.m'}, ...
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
    T = []; WV = [];
else
    [LoadingEngineDirectory, LoadingEngine, ext] = fileparts(LoadingEngine);
    MClust_NeuralLoadingFunction = LoadingEngine; 
    switch nargin
    case 1
        if nargout == 1
            T = feval(LoadingEngine, fn);
        else
            [T,WV] = feval(LoadingEngine, fn);
        end
    case 3
        [T,WV] = feval(LoadingEngine, fn, records_to_get, key);          
    otherwise
        error('Incorrect parameters passed to MClust_LoadNeuralData');
    end
end
    