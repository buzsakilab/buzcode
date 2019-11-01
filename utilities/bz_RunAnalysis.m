function [  ] = bz_RunAnalysis(analysisfunction,datasetfolder,varargin)
%bz_RunAnalysis(analysisfunction,datasetfolder) is a high level function that
%takes as its input an analysis function and returns the result of the
%analysis on multiple recordings
%
%INPUT
%   analysisfunction    A function name.  The function should be structured as
%                       [resultargs] = analysisfunction(basePath,figfolder)
%                       where resultargs are the variables from the
%                       analysis you'd like to save. (i.e. the results)
%   datasetfolder       The directory in which the dataset lives. Will find
%                       all basePaths in the directory.
%
%   (optional)
%       'savein'    'basePath' or 'functionPath' (default)
%       'cluster'   true/false, use true if running on the cluster to
%                   submit each recording as a batch job  NOT YET DONE
%       'basePath'  true if datasetfolder is a basePath, not a dataset folder
%
%
%TO ADD
%   option to select from spreadsheet of recordings with recording
%   properties
%   -option: save in one folder or save in the basePath
%   -option: fig folder
%
%DLevenstein 2017,2018
%%
p = inputParser;
addParameter(p,'savein','functionPath',@isstr)
addParameter(p,'cluster',false,@islogical)
addParameter(p,'basePath',false,@islogical)
parse(p,varargin{:})
savein = p.Results.savein;
cluster = p.Results.cluster;
BPin = p.Results.basePath;
%% Select Recordings to Analyze 

switch BPin
    case false
        [possiblebasePaths,possiblebaseNames] = bz_FindBasePaths(datasetfolder);
        [s,v] = listdlg('PromptString','Which recording(s) would you like analyze?',...
                        'ListString',possiblebaseNames);
        baseName = possiblebaseNames(s);
        basePath = possiblebasePaths(s); 
    case true
        basePath{1} = datasetfolder;
        baseName{1} = bz_BasenameFromBasepath(basePath{1});
end


%How many results are in the analysis function?
numresults = nargout(analysisfunction);

%% Make a figure folder
functionpath = which(analysisfunction);
functionpath = fileparts(functionpath);
figfolder = [fullfile(functionpath,'AnalysisFigs',analysisfunction),'/'];
if ~exist(figfolder,'dir')
    mkdir(figfolder)
end

%% Loop through the recordings and run the analysis

numrecs = length(baseName);
display(['Running Analysis on Recordings (',num2str(numrecs),')'])
for rr = 1:numrecs
    display(['Recording: ',num2str(rr),' of ',num2str(numrecs),' - ',baseName{rr}])
    
    [inputNames, outputNames] = get_arg_names(fullfile(functionpath,[analysisfunction,'.m']));
    inputNames = inputNames{1};
    if numresults>0
        outputNames = outputNames{1}; 
    end
    
    [results{1:numresults}] = feval(analysisfunction,basePath{rr},figfolder);
        
    for ss = 1:numresults
        eval([outputNames{ss} '= results{ss};']);
    end
    
    %Add: save date, baseName
    switch savein
        case 'basePath'
            savefolder = basePath{rr}; %or basePath...
        case 'functionPath'
            savefolder = figfolder; %or basePath...
        otherwise
            error('savein must be ''basePath'' or ''functionPath''')
    end
    savename = fullfile(savefolder,[baseName{rr},'.AnalysisResults.',analysisfunction,'.mat']);
    if numresults>0
        save(savename,outputNames{:});
    end
    
    close all
end

disp('ANALYSIS COMPLETE. HUZZAH!')

end

