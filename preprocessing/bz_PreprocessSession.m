function  bz_PreprocessSession(varargin)

%         bz_PreprocessSession(varargin)

%   Master fucntion to run the basic pre-processing pipeline for an
%   individual sessions. Is based on sessionsPipeline.m but in this case
%   works on an individual session basis no in a folfer with multiple ones.
% 

% INPUTS
%   <options>       optional list of property-value pairs (see table below)
%   expPath        - Basepath for experiment. It contains all session
%                       folders. If not provided, promt.
%   analogCh       - List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
%   forceSum       - Force make folder summary (overwrite, if necessary). Default false.
%   cleanArtifacts - Remove artifacts from dat file. By default, if there is analogEv in folder, is true.
%   stateScore     - Run automatic brain state detection with SleepScoreMaster. Default true.
%   spikeSort      - Run automatic spike sorting using Kilosort. Default true.
%   getPos         - get tracking positions. Default true. 
%   runAnalysis    - run summary analysis using AnalysisBatchScrip. Default false.
%   analysisList   - Logical array to indicate analysis to perform, according to the following list: 
%                           1. Spike-waveform, autocorrelogram and spatial position 
%                           2. Psth from analog-in inputs
%                           3. Slow-waves psth
%                           4. Ripples psh
%                           5. Theta, gamma and HFO profile. Theta and
%                               gamma mod
%                   Example: [1 1 1 1 1] or 'all' make all. Default 'all'    
%   pullData       - Path for raw data. Look for not analized session to copy to the main folder basepath. To do...
%
%  HISTORY: 
%     - Created based on sessionsPipeline: AntonioFR, 5/20

%  TO DO:
%   - Verify that data format and alysis output are compatible with CellExplorer
%   - Include Kilosort2 support
%   - Improve auto-clustering routine 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',pwd,@isdir); % by default, current folder
addParameter(p,'analogCh',[],@isnumeric);
addParameter(p,'forceSum',false,@islogical);
addParameter(p,'stateScore',true,@islogical);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'getPos',true,@islogical);
addParameter(p,'runAnalysis',false,@islogical);
addParameter(p,'analysisList','all');
addParameter(p,'cleanArtifacts',false,@islogical);
% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

expPath = p.Results.expPath;
analogCh = p.Results.analogCh;
forceSum = p.Results.forceSum;
stateScore = p.Results.stateScore;
spikeSort = p.Results.spikeSort;
getPos = p.Results.getPos;
runAnalysis = p.Results.runAnalysis;
analysisList = p.Results.analysisList;
cleanArtifacts = p.Results.cleanArtifacts;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
cd(expPath);

%% deals with xml.
if strcmp(expPath(end),filesep)
    expPath = expPath(1:end-1);
end
[~,basename] = fileparts(expPath);

disp('Check xml...');
if isempty(dir([basename '.xml'])) && isempty(dir('global.xml'))
    disp('No xml global file! Looking for it...');
    allpath = strsplit(genpath(expPath),';'); % all folders
    xmlFile = []; ii = 2;
    while isempty(xmlFile) && ii < size(allpath,2)
        disp(ii);
        cd(allpath{ii});
        xmlFile = dir('*.xml');
        ii = ii + 1;
    end
    if isempty(xmlFile)    
        [file, path] = uigetfile('*.xml','Select global xml file');
        copyfile(strcat(path,file),[basename '.xml']);
    else
        copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,basename,'.xml'));
    end
    cd(expPath);
end

%% Concatenate sessions
cd(expPath);
disp('Concatenate session folders...');
bz_ConcatenateDats(pwd,0,1);

%% Make SessionInfo
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  
save(strcat(basename,'.sessionInfo.mat'),'sessionInfo');
try
session = sessionTemplate(pwd,'showGUI',false); % look for rfh file - problem
save([basename '.session.mat'],'session');
catch
    warning('it seems that CellExplorer is not on your path');
end

%% Get analog and digital pulses
if  ~isempty(analogCh)
    [pulses] = bz_getAnalogPulses('analogCh',analogCh);
end
if ~exist('digitalIn.mat')
    digitalIn = bz_getDigitalIn('all','fs',sessionInfo.rates.wideband); 
end

%% Remove stimulation artifacts
if cleanArtifacts && ~isempty(analogCh)
    cleanPulses(pulses.ints{1}(:),'ch',0:session.extracellular.nChannels-mod(session.extracellular.nChannels,16)-1);
end

%% Make LFP
if isempty(dir('*.lfp'))
    try 
        bz_LFPfromDat(pwd,'outFs',1250); % generating lfp, NOTE: THIS FUNCTION WILL GENERATE A SESSIONINFO FILE!! WE NEED TO FIX THIS
        disp('Making LFP file ...');
    catch
        disp('Problems with bz_LFPfromDat, resampling...');
        ResampleBinary(strcat(sessionInfo.session.name,'.dat'),strcat(sessionInfo.session.name,'.lfp'),...
        sessionInfo.nChannels,1,sessionInfo.rates.wideband/sessionInfo.rates.lfp);
    end
end

%% Get brain states
if stateScore 
    try 
        if exist('pulses','var')
            SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
        else
            SleepScoreMaster(pwd,'noPrompts',true); % try to sleep score
        end
    catch
        disp('Problem with SleepScore skyping...');
    end
end

%% Kilosort concatenated sessions
if spikeSort
    disp('Spike sort concatenated sessions...');
        if  isempty(dir('*Kilosort*')) % if not kilosorted yet
        fprintf(' ** Kilosorting session %3.i of %3.i...');
            KiloSortWrapper;
            kilosortFolder = dir('*Kilosort*');
            try PhyAutoClustering(strcat(kilosortFolder.folder,'\',kilosortFolder.name)); % autoclustering
            catch
                warning('PhyAutoClustering not possible!!');
            end
            if exist('phyLink') && ~isempty(phyLink) % move phy link to
                kilosort_path = dir('*Kilosort*');
                try copyfile(phyLink, strcat(kilosort_path.name,filesep,'LaunchPhy')); % copy pulTime to kilosort folder
                end
            end
        end
    cell_metrics = ProcessCellMetrics('session', session);
end

%% Get tracking positions 
if getPos
    getSessionTracking;
end

%% deals with analysis list
if ischar(analysisList)
    if strcmpi(analysisList,'all')
        analysisList = [1 1 1 1 1 1 1];
        analysisList(2) = ~isempty(analogCh);
    else
        error('Analysis list format not recognized!');
    end
end

%% BatchAnalysis
if runAnalysis
    fprintf(' ** Summary %3.i of %3.i... \n',ii, size(allSess,1));
    if forceSum || (~isempty(dir('*Kilosort*')) && (isempty(dir('summ*')))) % is kilosorted but no summ
        disp('running summary analysis...');
        AnalysisBatchScript;         
    end
end

end

