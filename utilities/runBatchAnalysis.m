function [ results ] = runBatchAnalysis(analysisToRun,directoryList,fileList)
% runBatchAnalysis(analysisToRun,directoryList,fileList) is a high level function that
% takes as its input an analysis function, a list of folders, and a list of files.  It
% iterates through each folder, loading the listed files and running the desired analysis.
% It returns the results of the analysis in a struct with the same
% dimensions as directoryList.
%
%   INPUT
%   analysisToRun       A function name.  
%
%   inputArgs           A list of input arguments for the specified
%                       function to use
%                      
%   directoryList       A cell array that lists folders for which the 
%                       analysis should be run. 
%
%   fileList            A cell array that lists the files that must be
%                       loaded to run the function "analysisToRun
% 
%
%
% verbally designed by Dan Levenstein, Brendon Watson, David Tingley, 02/2017

% for functionpath and nargout of builtins we don't want the .m
if strcmp(analysisToRun(end-1:end),'.m')
    analysisToRun = analysisToRun(1:end-2);
end
%% Loop through the recordings and run the analysis
numresults = nargout(analysisToRun);
numrecs = length(directoryList);
functionpath = which([analysisToRun '.m']);
functionpath = fileparts(functionpath);


display(['Running Analysis: ' analysisToRun ', on ' num2str(numrecs) ' Recordings..'])
for rr = 1:numrecs
    display(['Recording: ',num2str(rr),' of ',num2str(numrecs),' - ',directoryList{rr}])
    try
        % load necessary files into workspace
        for i = 1:length(fileList)
            load([directoryList{rr} '/' fileList{i}])
        end
    catch 
       error(['the necessary *.mat file could not be found for ' directoryList{rr}]) 
    end
    
    [inputNames, outputNames] = get_arg_names(fullfile(functionpath,[analysisToRun '.m']));
    
    outputNames = outputNames{1}; inputNames = inputNames{1};
    for i=1:length(inputNames)
        inputVars{i} = eval(inputNames{i});
    end
    % generate function handle
    s = split(analysisToRun,'.');
    s = s{1};
    handle = eval(['@' s]);
    resTemp = feval(handle,inputVars{:});
    for i=1:length(outputNames)
       eval(['results(rr). ' outputNames{i} ' = resTemp'])
    end
end
end
