
function sessionsPipeline(varargin)

% expPipeline(varargin)

% Master script for launching the basic pre-processing of Intan recorded ephys data.
% Is meant to be run on a folder with multiple sesions from different
% animals and organizes them into same-day recordings as a first step. 
%%%%%%% EXPLAIN HERE NAME CONVENTION ASSUMED ...
% Also makes premilinary descriptive analysis and figures as an overview of the session.

%   1. Organizes data folders by sessions (sessions being all recordings from same day).
%   2. Concatenate sessions data.
%   3. Spike sort by kilosort.
%   4. Autocluster.
%   5. Makes a folder summary with spike-waveforms, autocorrelogram and spatial position; psth from analog-in inputs and
%       psth around slow waves and ripples. It requires AnalysisBatchScript.m

% INPUTS
%   <options>       optional list of property-value pairs (see table below)
%   expPath        - Basepath for experiment. It contains all session
%                       folders. If not provided, promt.
%   analogCh       - List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
%   forceSum       - Force make folder summary (overwrite, if necessary). Default false.
%   cleanArtifacts - Remove artifacts from dat file. By default, if there is analogEv in folder, is true.
%   spikeSort      - Run automatic spike sorting using Kilosort. Default true.  
%   pullData       - Path for raw data. Look for not analized session to copy to the main folder basepath. To do...
%
%  HISTORY: 
%   - Manu Valero-BuzsakiLab 2019
%   - Some reorganization, added state scoring and LFP creation, other minor changes: 5/20, AntonioFR
%
%  TO DO:
%   - Verify that data format and alysis output are compatible with CellExplorer
%   - Include Kilosort2 support
%   - Improve auto-clustering routine 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isdir);
addParameter(p,'analogCh',[],@isnumeric);
addParameter(p,'forceSum',false,@islogical);
addParameter(p,'stateScore',false,@islogical);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'analisysList','all');
addParameter(p,'cleanArtifacts',false,@islogical);
% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

expPath = p.Results.expPath;
analogCh = p.Results.analogCh;
forceSum = p.Results.forceSum;
stateScore = p.Results.stateScore;
spikeSort = p.Results.spikeSort;
analisysList = p.Results.analisysList;
cleanArtifacts = p.Results.cleanArtifacts;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});

%% deals with xml.
disp('Check xml...');
if isempty(dir('global.xml')) 
    disp('No xml global file! Looking for it...');
    xmlFile = []; ii = 2;
    while isempty(xmlFile) && ii < size(allpath,2)
        disp(ii);
        cd(allpath{ii});
        xmlFile = dir('*.xml');
        ii = ii + 1;
    end
    if isempty(xmlFile)    
        [file, path] = uigetfile('*.xml','Select global xml file');
        copyfile(strcat(path,file),'global.xml');
    else
        copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,'global.xml'));
    end
    cd(allpath{1});
end

%% Build sessions
disp('Building session folders (It asumes session as all folders recorded same day)...');
allFolder = dir(pwd);
for ii = 1:length(allFolder)
    if strlength(allFolder(ii).name) > 12 && isfolder(allFolder(ii).name) % if looks like a data folder
        folderFiels = strsplit(allFolder(ii).name,'_');
        if isempty(dir(strcat('*',folderFiels{2},'_','sess*'))) % if there is no session folder yet
            mkdir(strcat(folderFiels{1},'_',folderFiels{2},'_','sess',num2str(size(dir('*sess*'),1) + 1))); % create folder
        end
        if ~contains(folderFiels{3},'sess') % if it is not a session folder
            targetfoder = dir(strcat(folderFiels{1},'_',folderFiels{2},'_','sess*'));
            movefile(strcat('*',folderFiels{2},'_',folderFiels{3}),targetfoder.name); % move to session folder
        end
    end
end

% If recordered with intan Buz Edit, change dat file name
allpath = strsplit(genpath(expPath),';'); % all folders again
for ii = 1:size(allpath,2)
    if strlength(allpath{ii}) > 12 && isfolder(allpath{ii}) % if looks like a data folder
        cd(allpath{ii});
        if ~isempty(dir('amplifier_analogin_auxiliary_int16.dat')) % if recordered with intan Buz Edit
            movefile('amplifier_analogin_auxiliary_int16.dat','amplifier.dat');
            try  movefile('amplifier_analogin_auxiliary_int16.xml','amplifier.xml');
                movefile('amplifier_analogin_auxiliary_int16.nrs','amplifier.nrs');
            end
        end
    end
end

%% Concatenate sessions
cd(allpath{1});
allSess = dir('*_sess*');
disp('Concatenate session folders...');
for ii = 1:size(allSess,1)
    fprintf(' ** Concatenating session %3.i of %3.i... \n',ii, size(allSess,1));
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    delete(strcat(allSess(ii).name,'.xml'));% bring xml file
    copyfile(strcat(allpath{1},'\global.xml'),strcat(allSess(ii).name,'.xml'),'f');
    % [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    bz_ConcatenateDats(pwd,0,1); % concat files according to time of recording
end

%% Make SessionInfo
% [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  
% save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
session = sessionTemplate(pwd,'showGUI',false);

%% Remove stimulation artifacts
if cleanArtifacts && ~isempty(analogCh)
    [pulses] = bz_getAnalogPulses('analogCh',analogCh);
    cleanPulses(pulses.ints{1}(:));
end

%% Make LFP
if isempty(dir('*.lfp'))
    filename = split(pwd,'\'); filename = filename{end};
    ResampleBinary(strcat(filename,'.dat'),strcat(filename,'.lfp'),...
        session.extracellular.nChannels,1, session.extracellular.sr/session.extracellular.srLfp);
end

%% Get brain states
if stateScore %&& exist('StateScoreFigures','dir')~=7
    try 
        SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
    catch
        disp('Problem with SleepScore skyping...');
    end
end

%% Kilosort concatenated sessions
if spikeSort
    disp('Spike sort concatenated sessions...');
    for ii = 1:size(allSess,1)
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        if  isempty(dir('*Kilosort*')) % if not kilosorted yet
        fprintf(' ** Kilosorting session %3.i of %3.i... \n',ii, size(allSess,1));
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
    end
end

%% BatchAnalysis
for ii = 1:size(allSess,1)
    fprintf(' ** Summary %3.i of %3.i... \n',ii, size(allSess,1));
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    if forceSum || (~isempty(dir('*Kilosort*')) && (isempty(dir('summ*')))) % is kilosorted but no summ
        disp('running summary analysis...');
        sessionSummary;         
    end
end

end
