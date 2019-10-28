%function StateEditorProBeta(baseName, inputData, supressGUI, makePortable)
% SIMPLE USAGE (NO INPUT VARIABLES, if run from dir with .xml and .eeg/lfp):
% Add StateEditor to your Matlab path. Within Matlab navigate to a folder
% containing '.xml' and '.eeg' (or '.lfp') files. Type in the name of this 
% .m file and press enter (StateEditor will get the relevant baseName from 
% the name of the first available '.xml' file in the folder). The StateEditor
% loading GUI will guide you through channel selection and processing. Once 
% the TheStateEditor GUI loads, press 'H' for further information on using 
% StateEditor.
% 
% THE '.eegstates.mat' FILE: When it first runs on a new folder
% StateEditor creates a 'baseName.eegstates.mat' file. This file contains
% your channel selection as well as the (whitened) spectrograms. Subsequent
% runs on this folder will automatically load and use the selections found
% in the '.eegstates.mat' file, substantially speeding up the loading
% process. In order to view different channels you must first delete or
% rename this file. Note that in order to save space StateEditor does not
% save the lfp channels selected in the 'eegstates.mat' file unless the
% 'makePortable' input variable is set to 1.
% 
% SAVING AND LOADING STATE EDITOR WORK: Pressing 'S' allows you to save
% your StateEditor work. The default output file name is:
% 'baseName-states.mat'. By default this file contains a structure with 3
% fields: 
%    'states' - the state vector is of length N, where N is the number
%               of seconds	bins in your spectrogram 
%               (N = round((length(eeg)/eegFS) - 1). It has a value between
%               0 and 5 for each bin (0 = 'no state', 1 = 'awake', 2 = 
%               'Light/Drowzy', 3 = 'NREM', 4 = 'Intermediate', 5 = 'REM').
%   'events' -  a N by 2 matrix where the first column are event type ID's 
%               (1 through max of 10 types) and the 2nd column are the 
%               event times in seconds (in arbitrary precision). 
%   'transitions' - a Nx3 matrix of exact state transition times. 1st 
%                   collumn: state number, 2nd collumn: state start time 
%                   (in seconds), 3rd column: 	end time (in seconds). The 
%                   transition matrix is a higher resolution complement to
%                   the state vector for those who wish to choose their
%                   states with a resolution greater than 1Hz. In order to 
%                   load states or events into StateEditor save a structure 
%                   formatted as specified above.
% 
% LOADING VARIABLES DIRECTLY FROM MATLAB: The structure 'inputData' allows
% for loading directly from existing matlab variables (NOTE: the details
% below are only relevant to those that do not have access to .lfp/.eeg
% files or wish to bypass the StateEditor loading GUI) It must contain the
% following (case-sensitive) fields: 
%       'rawEeg': this is a cell array of size N where N is the number of 
%               eeg channels to load (maximum of 3). Each cell of 'rawEeg' 
%               must be a vector 1 eeg channel of length n where n is the
%               number of samples 'Chs': 1xN matrix of channel numbers 
%               where N is the number of lfp channels to load 
%               (for instance inputData.Chs = [20, 39])
%       'eegFS': lfp sampling frequency (default: 1250 Hz)
%       'Chs': List of channel numbers corresponding to rawEeg
%       'nChs': number of input channels in rawEeeg
%       'MotionType':
%           Must be one of the following (case-sensitive) strings:
%           -'none'
%           -'Whl' 
%           -'Channels (accelerometer)' 
%           -'Channels (MEG)' 
%           -'File'
%       'motion': a 1xN vector of movement data to be displayed in the 
%           motion panel. N is the number of second bins in the spectrogram
%           of the lfp channels and has the value: 
%               round((length(rawEeg{1}/eegFS) - 1). 
%   NOTE: to use a pre-processed motion signal, set the 'MotionType' to 'none' and
% pass the actual motion signal in the field 'motion' as described above.
% Conversely, if you wish to pass a motion signal to be processed (for
% instance, accelerometer or MEG channel(s), pass these channels in through
% the field 'MotionSignal' (which must be yxn vector or matrix (where y is
% the number of motionsignal channels, and n is the number of samples in
% eegFS samples/second)) and set the 'MotionType' to the desired processing
% type.
% 
% 
% SUPRESSING THE GUI LOADER AND THE STATE EDITOR GUI: 
% If you wish to bypass the StateEditor GUI loader use the 'inputData' 
% structure to pass StateEditor the variables: 'Chs' and 'MotionType'. If 
% applicable also  select the motion signal channels through the variable 
% 'mChs'. Note that StateEditor must still be called from the folder in 
% which the relevant '.xml' and '.eeg/.lfp' files are stored. In order to 
% suppress the StateEditor GUI, set the supressGUI variable to 1. This can 
% be useful for those wishing to pre-create the '.eegstates.mat' file as 
% part of their data pre-processing.
%
% Dependency: Does require LoadXml.m from XmlTree
% http://www.artefact.tk/software/matlab/xml/
%
%created by Andres Grosmark at Gyuri Buzsaki's lab, 12/2012.
%Improvements by Brendon Watson, DLevenstein
%Many subfunctions, mostly from Anton Sirota, but also from Adrien Peyrache
%and others have been included as subfunctions of this script to reduce
%dependency issues.
%



function TheStateEditor(baseName, inputData, supressGUI, makePortable)

%% get baseName if doesn't exist, save

if exist('inputData', 'var');
    if ~isempty(inputData);
        fields = fieldnames(inputData);
        for i = 1:length(fields)
            eval([fields{i} ' = inputData.', fields{i}, ';']);
        end
    end
end

if exist('Chs', 'var') & exist('rawEeg', 'var') & exist('MotionType', 'var') & ~exist('nCh', 'var')
    nCh = max(Chs);
end

if ~exist('eegFS', 'var')
        eegFS = 1250;%this is dangerous, creates some problems I'll try to fix below, BW
end
LoadFromPortable = 0;
if ~exist('baseName','var')
    baseName = [];
end
if isempty(baseName)
    xmlBase = dir('*xml');
    if length(xmlBase) == 0
%         if FileExistsIn('*.eegstates.mat')
%             d = dir('*.eegstates.mat');
%             StateInfo = load(d(1).name);
%             StateInfo = StateInfo.StateInfo;
%             
%             if ~isfield(StateInfo, 'rawEeg')
%                 warndlg({['Only ''.eegstates.mat'' file containing '],
%                     ['the eeg/lfp signals can be loaded without'],
%                     ['an ''.xml'' file present.']});
%                 return;
%             else
%                 d = dir('*eegstates.mat');
%                 baseName = d(1).name(1:(end - 14));
%                 LoadFromPortable = 1;
%             end
%         else
%             
%             warndlg('No ''*xml'' file found. Quitting now. Bye bye.');
%             return
%             
%         end
    else
        xmlBase = xmlBase(1);
        choice = questdlg(['No basename entered, use ', xmlBase.name(1:(end - 4)), ' as file basename?'],'No Basename','Yes','Cancel','Yes');
        if strmatch(choice,'Cancel')
            return
        elseif strmatch(choice,'Yes')
            baseName = xmlBase.name(1:(end - 4));
        end
    end
end


if FileExistsIn([baseName, '.eegstates.mat']);
    StateInfo = load([baseName, '.eegstates.mat']);
    StateInfo = StateInfo.StateInfo;
    if isfield(StateInfo, 'rawEeg')
        LoadFromPortable = 1;
    end
end
suffix = [];


%         if FileExistsIn('*.eegstates.mat')
%             d = dir('*.eegstates.mat');
%             StateInfo = load(d(1).name);
%             StateInfo = StateInfo.StateInfo;
%             
%             if ~isfield(StateInfo, 'rawEeg')
%                 warndlg({['Only ''.eegstates.mat'' file containing '],
%                     ['the eeg/lfp signals can be loaded without'],
%                     ['an ''.xml'' file present.']});
%                 return;
%             else
%                 d = dir('*eegstates.mat');
%                 baseName = d(1).name(1:(end - 14));
%                 LoadFromPortable = 1;
%             end
%         else
%             
%             warndlg('No ''*xml'' file found. Quitting now. Bye bye.');
%             return
%             
%         end







if ~exist('supressGUI', 'var')
    supressGUI = 0;
end

if exist('Chs', 'var') & exist('MotionType', 'var')
    supressLoadGUI = 1;
else
    supressLoadGUI = 0;
end

if ~exist('makePortable', 'var')
    makePortable = 0;
end

if ~exist('spikeInfo', 'var')
    spikeInfo = 0;
end


%These parameters are passed through all functions
FO.downsampleGoal = 312.5;% display Hz goal, to save memory... will calculate downsample factor to match (ie 4 if 1250hz lfp file)
FO.eegAlreadyDownsampled = 0; %will have been downsampled if saved by user as raw
FO.baseName = baseName;  %Includes basePath
FO.basePath = fileparts(baseName);
if isempty(FO.basePath)
    FO.basePath = pwd;     %basePath is assumed to be pwd... 
end
FO.eegDisplaySeconds = 2; %show 2 seconds of eeg
FO.maxFreq = 40; %default starting frequency extent
FO.hanningW = 10; %default hanning smoothing window
FO.EegUpdateBoolean = logical(0);
FO.stateAxisToggleModeBool = logical(1);

%FO.lax - handle for state label axes
%FO.ilab - handle for state label plot
%FO.sax{1:nCh} - handles for the spectrogram axes
%FO.iSpec{1:nCh} - handles for the spectrogram objects
%FO.max - handle for motion axes
%FO.Mplot - handle for motion plot
%FO.eax{1:nCh} - handles for eeg axes
%FO.Eplot{1:nCh} - handles for eeg plots
%FO.sMiddline{1:nCh} - handles for spectrogram center of view line objects
%FO.mMidline - handle for motion center of view line object
%FO.eMidline{1:nCh} - handle for eeg center of view line objects


%% account for possibility that some lft files are .eeg and some are .lft
if ~(exist('rawEeg', 'var') & exist('Chs', 'var') & exist('nCh', 'var') & exist('MotionType', 'var'))
    if LoadFromPortable == 0
        if FileExistsIn([baseName,'.eeg'])
            suffix = '.eeg';
        else
            if FileExistsIn([baseName,'.lfp'])
                suffix = '.lfp';
            else
% %                 try 
% %                     basepath = cd;
% %                     eeglfppath = findsessioneeglfpfile(baseName,basepath);
% %                 catch
%                     disp(['Error: ', baseName, '.eeg or .lfp not found.'])
%                     disp(['Quitting now. Bye bye.']);
%                     return
% %                 end
            end
        end
        
        if ~FileExistsIn([baseName, '.xml'])
            disp(['Error: ', baseName, '.xml not found.'])
            disp(['Quitting now. Bye bye.']);
            return;
        end
    end
end


%% check for prior processing
states = [];
if FileExistsIn([baseName,'.eegstates.mat'])
    if ~exist('StateInfo','var')
        StateInfo = load([baseName,'.eegstates.mat']);
        StateInfo = StateInfo.StateInfo;
    end
    eegloaded = 0;

    if isfield(StateInfo, 'rawEeg')%if made portable, take from StateInfo
        rawEeg = StateInfo.rawEeg;
        eegloaded = 1;
    elseif FileExistsIn([baseName,'.RawEEG.eegstates.mat'])%if raw eeg saved by user (separate file)
        load(fullfile(FO.basePath,[FO.baseName, '.RawEEG.eegstates.mat']))
        %check if channels are right
        schan = StateInfo.Chs(:);
        echan = RawEegData.channels(:);
        if length(echan) == length(schan)
            if sum(abs(echan-schan)) == 0
                rawEeg = RawEegData.data;
                for eix = 1:length(rawEeg)
                    rawEeg{eix} = double(rawEeg{eix});
                end
                eegloaded = 1;
                FO.eegAlreadyDownsampled = 1;
            end
        end
        clear RawEegData schan echan
    end
    
    if ~eegloaded %load from .lfp/eeg if not loaded in above if/elseif
        if ~exist('rawEeg', 'var')
            rawEeg = {};
            if isfield(StateInfo,'eegFS')
                eegFS = StateInfo.eegFS;%this was missing and caused probs with default 1250Hz assumption if data not at 1250hz
            else %allows compatibility with old files...
                eegFS = 1250;
                StateInfo.eegFS = eegFS;%...should fix them too
            end
            Chs = StateInfo.Chs;
            nCh = StateInfo.nCh;
            disp([baseName, '.eegstates.mat loaded']);
            disp(['Using channel(s): ', int2str(Chs)]);
            
            disp('Retrieving eeg channel(s)...');
            rawEeg = {};
            eeg = [];
            try
                for i = 1:length(Chs)
                    e = LoadChanArch(baseName, Chs(i));
                    eeg = [eeg; e];
                    e = [];
                end
            catch
                disp('No eeg in your eegstates.mat. Loading from .lfp/.eeg file...');
                try %first try bz_getLFP
                    eeg = bz_GetLFP(Chs,'basepath',FO.basePath,'noPrompts',true);
                    eeg = single(eeg.data)';
                catch
                    try
                        % try Anton's LoadBinary
                        eeg = LoadBinary([baseName, suffix], Chs+1, nCh, [], 'int16', 'single');
                    catch
                        %Otherwise try to use Micheal Zugaro
                        eeg = LoadBinaryIn([baseName, suffix], 'channels', Chs+1, 'nChannels', nCh)';
                        eeg = single(eeg);
                    end
                end
                
            end
            
            for i = 1:length(Chs)
                rawEeg{i} = eeg(i, :);
            end
            disp('Done.');
        end
    end

else
    StateInfo = [];
    info1 = bz_getSessionInfo(FO.basePath,'noPrompts',true);
    eegFS = info1.lfpSampleRate;

    if ~exist('nCh', 'var')
%             info1 = LoadXmlIn([baseName, '.xml']);
            nCh = info1.nChannels;
    end
    if supressLoadGUI == 0
        
        global answer1
        answer1 = 0;
        
        if exist([baseName '.SleepScoreLFP.LFP.mat'],'file')
            load([baseName '.SleepScoreLFP.LFP.mat'])
            defaultchans = ([SleepScoreLFP.SWchanID SleepScoreLFP.THchanID]);
            if length(defaultchans)<2
                defaultchans = num2str(defaultchans);
            elseif length(defaultchans)>=2
                dcout = num2str(defaultchans(1));
                for chanidx = 2
                    dcout = strcat(dcout,',',num2str(defaultchans(chanidx)));
                end
                defaultchans = dcout;
            end
        else
            defaultchans = '';
        end
        
        inputFig = figure('Position', [280   453   550   250], 'MenuBar', 'none', 'numbertitle', 'off', 'name', [baseName, ' channel selection']);
        warning('off', 'MATLAB:hg:default_child_strategy:IllegalPermutation')
        annotation('textbox',  'Position', [0.02, 0.87, 0.9, 0.07], 'string', ['\bf\fontsize{10}Please choose up to 3 eeg channels (base 0, nCh = ', int2str(nCh), '):'], 'EdgeColor', 'none');
        Channels1 = uicontrol('style', 'edit', 'FontSize', 10, 'Units', 'Normalized', 'Position', [0.37, 0.75, 0.45, 0.08],'String',defaultchans);
        set(Channels1, 'HorizontalAlignment', 'left');
        
        
        annotation('textbox',  'Position', [0.02, 0.65, 0.9, 0.07], 'string', ['\bf\fontsize{10}Choose a motion signal to use:'], 'EdgeColor', 'none');
        mOptions = {'.EMGFromLFP.LFP.mat';'Load From .whl file (head tracking)';'Load from eeg ch(s) (accelerometer/motion pad)';'Load from eeg ch(s) (MEG)';'Load from .mat file vector';'Load from TimeValue Pair .mat';'None'};
        mInput = uicontrol('style', 'popupmenu', 'string', mOptions );
        set(mInput, 'Units', 'normalized', 'Position', [0.37, 0.5, 0.62, 0.1]);
        
        annotation('textbox',  'Position', [0.02, 0.405, 0.9, 0.07], 'string', ['\bf\fontsize{10}Choose motion signal channel(s) (if applicable):'], 'EdgeColor', 'none');
        
        mChannels1 = uicontrol('style', 'edit', 'FontSize', 10, 'Units', 'Normalized', 'Position', [0.37, 0.295, 0.45, 0.08]);
        set(mChannels1, 'HorizontalAlignment', 'left');
        
        proceed1 = uicontrol('style', 'pushbutton', 'string', 'Go!', 'Callback', 'global answer1; uiresume(gcbf); answer1 = 1;');
        set(proceed1, 'Units', 'normalized', 'Position', [0.4, 0.05, 0.25, 0.12], 'FontSize', 12);
        
        cancel1 = uicontrol('style', 'pushbutton', 'string', 'Cancel', 'Callback', 'global answer1; uiresume(gcbf); answer1 = 0;');
        set(cancel1, 'Units', 'normalized', 'Position', [0.7, 0.05, 0.25, 0.12], 'FontSize', 12);
        
        uiwait(inputFig);
        Chs = get(Channels1, 'string');
        mIn = mOptions{get(mInput, 'Value')};
        mChs = get(mChannels1, 'string');
        
       
        clf(inputFig);
        close(inputFig);
        if answer1 == 0;
            return;
        end
    end
    if ischar(Chs)
        Chs = str2num(Chs);
    end
%     if exist('MotionType', 'var') & ~exist('mIn', 'var')
%         switch MotionType
%             case 'none'
%                 mIn = 1;
%             case  'Whl'
%                 mIn = 2;
%             case 'Channels (accelerometer)'
%                 mIn = 3;
%             case 'Channels (MEG)'
%                 mIn = 4';
%             case 'File'
%                 mIn = 5;
%         end
%     end
    
    
    if ~iscell(Chs)
        if sum(Chs >= 0 & Chs <= nCh) ~= length(Chs) | isempty(Chs)
            b = msgbox('Error: Incorrect channel selection. Quiting now. Bye bye.');
            uiwait(b);
            return;
        end
    end
    
    
    if ~exist('rawEeg', 'var')
        weeg = {};
        fspec = {};
        
        disp(['Loading eeg channels: ', int2str(Chs)]);
        try %first try bz_getLFP
            eeg1 = bz_GetLFP(Chs,'basepath',FO.basePath,'noPrompts',true);
            eeg1 = single(eeg1.data)';
        catch
            try
                % try Anton's LoadBinary
                eeg1 = LoadBinary([baseName, suffix], Chs+1, nCh, [], 'int16', 'single');
            catch
                %Otherwise try to use Micheal Zugaro
                eeg1 = LoadBinaryIn([baseName, suffix], 'channels', Chs+1, 'nChannels', nCh)';
                eeg1 = single(eeg1);
            end
        end
        disp('Done.');
        for i = 1:length(Chs)
            
            rawEeg{i} =eeg1(i, :);
            if iscell(Chs)
                disp(['Whitening and computing spectrogram for channel ', Chs{i},'. This will all be over in a minute.']);
            else
                disp(['Whitening and computing spectrogram for channel ', int2str(Chs(i)),'. This will all be over in a minute.']);
            end
            try
                fspec{i} = LoadSpecArch(baseName, [], Chs(i), 1, 0, 3072, [0 200], 1, []);
            catch
                fspec{i} =[];
                
                weeg{i} =  WhitenSignalIn(rawEeg{i},eegFS*2000,1);
                [fspec{i}.spec, fspec{i}.fo, fspec{i}.to] = mtchglongIn(weeg{i}, 3072, eegFS, eegFS, 0, [], [], [], [0 200]);
                fspec{i}.spec = single(fspec{i}.spec);
                fspec{i}.info.Ch = Chs(i);
                fspec{i}.info.FileInfo.name = [baseName, suffix];
            end
            disp('Done.');
            
        end
    else
        for i = 1:length(Chs)
            
            
            if iscell(Chs)
                disp(['Whitening and computing spectrogram for channel ', Chs{i}, '. This will all be over in a minute.']);
            else
                 disp(['Whitening and computing spectrogram for channel ', int2str(Chs(i)),'. This will all be over in a minute.']);
            end
            try
                fspec{i} = LoadSpecArch(baseName, [], Chs(i), 1, 0, 3072, [0 200], 1, []);
            catch
                fspec{i} =[];
                weeg{i} =  WhitenSignalIn(rawEeg{i}, eegFS*2000,1);
                [fspec{i}.spec, fspec{i}.fo, fspec{i}.to] = mtchglongIn(weeg{i}, 3072, eegFS, eegFS, 0, [], [], [], [0 200]);
                fspec{i}.spec = single(fspec{i}.spec);
                fspec{i}.info.Ch = Chs(i);
                fspec{i}.info.FileInfo.name = [baseName, suffix];
            end
            disp('Done.');
            
        end
        
    end
    
    
    if ~exist('mChs', 'var')
        mChs = [];
    end
    
    if ~exist('motion', 'var')
        if ~isempty(mChs) & ischar(mChs)
            mChs = str2num(mChs);
        end
        
        
        if strcmp(mIn,'Load from eeg ch(s) (accelerometer/motion pad)')
        %         if mIn == 3
            if ischar(mChs)
                mChs  = str2num(mChs);
            end
            if (sum(mChs > 0 & mChs <= nCh) ~= length(mChs)) | isempty(mChs)
                b = msgbox('Error: Incorrect motion channel selection. Quiting now. Bye bye.');
                uiwait(b);
                return;
            end
        end
       
        
        MotionType = [];
        switch(mIn)
            case '.EMGFromLFP.LFP.mat'
                motion = [];
                MotionType = 'EMGFromLFP';
                if exist([baseName,'.EMGFromLFP.LFP.mat'],'file')
                    tpath = [baseName,'.EMGFromLFP.LFP.mat'];
                else
                    [name, path] = uigetfile('*.mat', 'EMG: Choose a file with time:val pairs:');
                    tpath = fullfile(path,name);
                end
                load(tpath)%should now have the EMG variable with fields 
                
                tos = fspec{1}.to;
                times = EMGFromLFP.timestamps;
                vals = EMGFromLFP.data;
                motion = ResampleTolerant_IN(vals,length(tos),length(times));
            case 'None'
                motion = [];
                MotionType = 'none';
            case 'Load From .whl file (head tracking)'
                disp('Loading and preprocessing motion data from .whl file...');
                motion = LoadFromWhl(baseName, fspec{1}.to);
                if sum(isnan(motion)) ~= 0
                    disp(['Note that ', num2str(mean(isnan(motion))*100), '% of the motion values are NaNs']);
                    disp('Proceeding...');
                end
                MotionType = 'Whl';
                mChs = [];
                disp('Done.');
            case 'Load from eeg ch(s) (accelerometer/motion pad)'
                disp(['Loading and preprocessing motion data from channel(s) ', int2str(mChs), '...']);
                if exist('motionSignal', 'var')
                    meeg = motionSignal;
                else
                    try
                        %First try Anton's LoadBinary
                        meeg = LoadBinary([baseName, suffix], mChs+1, nCh, [], 'int16', 'single');
                    catch
                        %Otherwise try to use Micheal Zugaro
                        meeg = LoadBinaryIn([baseName, suffix], 'channels', mChs+1, 'nChannels', nCh)';
                        meeg = single(meeg);
                    end
                end
                meeg = abs(zscore(meeg')');
                meeg = sum(meeg, 1);
                forder = 500;
                forder = ceil(forder/2)*2;
                EEGSR = eegFS;
                lowband = 0.1;
                highband = 1;
                firfiltb = fir1(forder,[lowband/EEGSR*2,highband/EEGSR*2]);
                meeg = filter2(firfiltb,  meeg);
                motion = mean(reshape(meeg(1:(length(meeg) - mod(length(meeg), eegFS))), eegFS, []), 1);
                if length(motion) == (length(fspec{1}.to) + 1)
                    motion = motion(1:(end - 1));
                end
                MotionType = 'Channels (accelerometer)';
                disp('Done.');
            case 'Load from eeg ch(s) (MEG)'
                disp(['Loading and preprocessing meg data from channel(s) ', int2str(mChs), '...']);
                if exist('motionSignal', 'var')
                    meeg = motionSignal;
                else
                    try
                        %First try Anton's LoadBinary
                        meeg = LoadBinary([baseName, suffix], mChs+1, nCh, [], 'int16', 'single');
                    catch
                        %Otherwise try to use Micheal Zugaro
                        meeg = LoadBinaryIn([baseName, suffix], 'channels', mChs+1, 'nChannels', nCh)';
                        meeg = single(meeg);
                    end
                end
                meeg = zscore(meeg')';
                meeg = sum(meeg, 1);
                forder = 500;
                forder = ceil(forder/2)*2;
                EEGSR = eegFS;
                lowband = 100;
                highband = 600;
                firfiltb = fir1(forder,[lowband/EEGSR*2,highband/EEGSR*2]);
                meeg = filter2(firfiltb,  meeg);
                meeg = zscore(meeg).^2;
                lowband = 0.1;
                highband = 1;
                firfiltb = fir1(forder,[lowband/EEGSR*2,highband/EEGSR*2]);
                meeg = filter2(firfiltb,  meeg);
                motion = mean(reshape(meeg(1:(length(meeg) - mod(length(meeg), eegFS))), eegFS, []), 1);
                if length(motion) == (length(fspec{1}.to) + 1)
                    motion = motion(1:(end - 1));
                end
                MotionType = 'Channels (MEG)';
                disp('Done.');
                
            case 'Load from TimeValue Pair .mat'
                MotionType = 'TimeVal';
                varname = [];
%                 b = msgbox(['Note: Motion vector must be 1xn where n is the number of time bins in seconds (n = ', int2str(length(fspec{1}.to)), ')']);
%                 uiwait(b);
                [name, path] = uigetfile('*.mat', 'Choose a file with time:val pairs:');
                matobj = matfile(fullfile(path,name));
                w = whos(matobj);
                if length(w)>1
                    for a = 1:length(w)
                        n{a} = w(a).name;
                    end
                    varname = listdlg('ListString',n,'SelectionMode','Single','Name','Variable choice','PromptString','Choose variable to load');
                    varname = n{varname};
                end
                
                tos = fspec{1}.to;
                if isempty(varname);
                    motion = LoadTimeStampValuePairs(tos,fullfile(path,name));
                else
                    motion = LoadTimeStampValuePairs(tos,fullfile(path,name),varname);
                end
                
            case 'Load from .mat file vector'
                MotionType = 'File';
                b = msgbox(['Note: Motion vector must be 1xn where n is the number of time bins in seconds (n = ', int2str(length(fspec{1}.to)), ')']);
                uiwait(b);
                [name, path] = uigetfile('*mat', 'Choose a motion vector to load:');
                motion = load([path, name]);
                if isstruct(motion)
                    f1 = fieldnames(motion);
                    motion = motion.(f1{1});
                end                    
                
                mChs = name;
        end
        
    end
    if ~exist('nCh', 'var')
        nCh = length(Chs);
    end
    StateInfo.nCh = nCh;
    StateInfo.Chs = Chs;
    StateInfo.mChs = mChs;
    StateInfo.MotionType = MotionType;
    StateInfo.fspec = fspec;
    StateInfo.motion = motion;
    StateInfo.eegFS = eegFS;
%     if makePortable == 1
%         StateInfo.rawEeg = rawEeg;
%     end
    
    disp(['Saving ', baseName, '.eegstates.mat...']);
    try
        save([baseName, '.eegstates.mat'], 'StateInfo');
    catch
        warndlg(['Failed to save ' , baseName, '.eegstates.mat']);
    end
end

%Load Previous state tagging in SleepState.states.mat
thispath = FO.basePath;
timevector = StateInfo.fspec{1}.to;
states = bz_LoadStates_StateEditorWrapper_In(thispath,timevector);

if eegFS>FO.downsampleGoal
    FO.downsample = round(eegFS/FO.downsampleGoal);
else
    FO.downsample = 1;
end

disp('So far so good. Now, loading StateEditor GUI. This is going to be great!');

if supressGUI == 1
    return;
else
    StateEditorSetup(StateInfo.fspec, StateInfo.motion, states, rawEeg, FO, eegFS);
end
% 
% if exist([baseName,'-states.mat'],'file')
%     LoadStatesAutoNoMsgs
% end


end


function StateEditorSetup(f, MP, States, rawEeg, FO, eegFS)

if ~iscell(f)
    a = f; e = rawEeg;
    f = {}; rawEeg = {};
    f{1} = a; rawEeg{1} = e;
    a = []; e =[];
end
nCh = length(f);
FO.nCh = nCh;

Chs =[];
for i = 1:length(f)
    Chs = [Chs; f{i}.info.Ch];
end
FO.Chs = Chs;

if isempty(States)
    States = zeros(1, length(f{1}.to));
end

for i = 1:nCh
    if FO.eegAlreadyDownsampled
        FO.eeg{i} = rawEeg{i};
    else
        FO.eeg{i} = rawEeg{i}(1:FO.downsample:end);
    %     FO.eeg{i} = (eeg{i}(1:FO.downsample:end)/2150)/1000;%why was this division?? To convert to volts in some old system?  
    end
end

FO.clickPoint = [];
FO.startLine = {};
FO.States = States;
FO.stateHistory = {};
FO.stateHistoryNum = 0;
FO.newStates = {};
FO.currAction = 'Browse';
FO.eegFS = eegFS;


FO.overlayLines = {};

posvar = get(0,'Screensize');
posvar(1) = 20;
posvar(2) = 20;
posvar(3) = posvar(3)-100;
posvar(4) = posvar(4)-100;
FO.fig = figure('KeyReleaseFcn', {@DefKey}, 'Position', posvar);
set(FO.fig, 'numbertitle', 'off', 'name', ['States: ', FO.baseName]);
set(FO.fig,'WindowButtonDownFcn', {@MouseClick}, 'WindowButtonUpFcn', {@unMouseClick}, 'Units', 'normalized');
set(FO.fig, 'WindowButtonMotionFcn', {@Nothing}, 'WindowScrollWheelFcn', {@MouseScroll});
set(FO.fig, 'CloseRequestFcn', {@CloseDialog});
set(FO.fig, 'Tag', 'StateEditorMaster');
FO.madeChanges = 0;
FO.startLocation = [];

%%%%%%%%%%%%%%Axes Positions%%%%%%%%%%%%%%%%%%
positions = [];
switch nCh
    case 1
        position.lax = [0.0500    0.940    0.8000    0.0500];
        position.sax{1} = [0.0500    0.3100    0.8000    0.6200];
        position.MP = [0.0500    0.170    0.8000    0.1000];
        position.eax{1} = [0.0500    0.04    0.8000    0.1000];
        position.eegCh{1} =   [0.8100    0.0350    0.1000    0.1000];
        
        a = annotation('textbox', 'Units', 'Normalized', 'EdgeColor', 'none');
        if iscell(FO.Chs)
            set(a, 'String', ['\bf\color{black}\fontsize{11}Ch ', FO.Chs{1}], 'Position', [-0.0050    0.7800    0.1000    0.1000]);
        else
            set(a, 'String', ['\bf\color{black}\fontsize{11}Ch', int2str(FO.Chs(1))], 'Position', [-0.0050    0.7800    0.1000    0.1000]);
        end
%         annotation('textarrow',[0.02 0.02],[0.14 0.14],'string','\fontsize{10} EEG (m.V.)', ...
%             'HeadStyle','none','LineStyle', 'none', 'TextRotation',90);
        
        position.eegWidth = [0.05, 0.0350, 0.1, 0.1];
    case 2
        position.lax = [0.0500    0.940    0.8000    0.0500];
        position.sax{1} = [0.05         0.62          0.8         0.31];
        position.sax{2} = [0.05         0.3          0.8         0.31];
        position.MP = [0.05         0.195          0.8          0.07];
        position.eax{1} =  [0.05         0.105         0.8          0.06];
        position.eax{2} =  [0.05         0.04          0.8          0.06];
        position.eegCh{1} =    [0.81, 0.069, 0.1, 0.1];
        position.eegCh{2} = [0.81, 0.004, 0.1, 0.1];
        position.eegWidth = [0.05, 0.069, 0.1, 0.1];
        a = annotation('textbox', 'Units', 'Normalized', 'EdgeColor', 'none');
        if iscell(FO.Chs)
            set(a, 'String', ['\bf\color{black}\fontsize{11}Ch ', FO.Chs{1}], 'Position', [-0.005, 0.79, 0.1, 0.1]);
        else
            set(a, 'String', ['\bf\color{black}\fontsize{11}Ch', int2str(FO.Chs(1))], 'Position', [-0.005, 0.79, 0.1, 0.1]);
        end
        
        a = annotation('textbox', 'Units', 'Normalized', 'EdgeColor', 'none');
        
        if iscell(FO.Chs)
            set(a, 'String', ['\bf\color{black}\fontsize{11}Ch', FO.Chs{2}], 'Position', [-0.005, 0.49, 0.1, 0.1]);
        else
            set(a, 'String', ['\bf\color{black}\fontsize{11}Ch', int2str(FO.Chs(2))], 'Position', [-0.005, 0.49, 0.1, 0.1]);
        end
        annotation('textarrow',[0.02 0.02],[0.14 0.14],'string','\fontsize{10} EEG (m.V.)', ...
            'HeadStyle','none','LineStyle', 'none', 'TextRotation',90);
    case 3
        position.lax = [0.0500    0.945    0.8000    0.0475];
        position.sax{1} = [0.0500    0.730    0.8000    0.21];
        position.sax{2} = [0.0500    0.515    0.8000    0.21];
        position.sax{3} = [0.0500    0.3100    0.8000    0.20];
        position.MP = [0.05        0.23          0.8         0.07];
        position.eax{1} = [0.05         0.15          0.8          0.05];
        position.eax{2} = [0.05         0.095          0.8          0.05];
        position.eax{3} = [0.05         0.04          0.8          0.05];
        position.eegWidth = [0.05, 0.105, 0.1, 0.1];
        
        position.eegCh{1} =   [0.8100    0.095    0.1000    0.1000];
        position.eegCh{2} = [0.8100    0.0400    0.100    0.1000];
        position.eegCh{3} = [0.8100    -0.015    0.1000    0.1000];
        a1 = annotation('textbox', 'Units', 'Normalized', 'EdgeColor', 'none');
        a2 = annotation('textbox', 'Units', 'Normalized', 'EdgeColor', 'none');
        a3 = annotation('textbox', 'Units', 'Normalized', 'EdgeColor', 'none');
        
        if iscell(FO.Chs)
            set(a1, 'String', ['\bf\color{black}\fontsize{11}Ch',FO.Chs{1}], 'Position', [-0.0050    0.8200    0.1000    0.1000]);
            set(a2, 'String', ['\bf\color{black}\fontsize{11}Ch', FO.Chs{2}], 'Position', [-0.0050    0.6150    0.1000    0.1000]);
            set(a3, 'String', ['\bf\color{black}\fontsize{11}Ch', FO.Chs{3}], 'Position', [-0.0050    0.410    0.1000    0.1000]);
        else
            set(a1, 'String', ['\bf\color{black}\fontsize{11}Ch', int2str(FO.Chs(1))], 'Position', [-0.0050    0.8200    0.1000    0.1000]);
            set(a2, 'String', ['\bf\color{black}\fontsize{11}Ch', int2str(FO.Chs(2))], 'Position', [-0.0050    0.6150    0.1000    0.1000]);
            set(a3, 'String', ['\bf\color{black}\fontsize{11}Ch', int2str(FO.Chs(3))], 'Position', [-0.0050    0.410    0.1000    0.1000]);
        end
        %         annotation('textarrow',[0.02 0.02],[0.14 0.14],'string','\fontsize{10} EEG (m.V.)', ...
        %             'HeadStyle','none','LineStyle', 'none', 'TextRotation',90);
end
FO.xplotLims = [0.0500, 0.85];

yplotLims = [];
yplotLFPLims = [];
for i = 1:FO.nCh
    yplotLims = [yplotLims; position.sax{i}(2), position.sax{i}(2) + position.sax{i}(4)];
    yplotLFPLims = [yplotLFPLims; position.eax{i}(2), position.eax{i}(2) + position.eegCh{i}(4)];
    
end
yplotLims = [yplotLims; position.lax(2), position.lax(2) + position.lax(4)];
yplotLims = [yplotLims; position.MP(2), position.MP(2) + position.MP(4)];
FO.yplotLims = yplotLims;
FO.yplotLFPLims = [yplotLFPLims];

colors = {};
FO.resolution = 0.5;
f1 = round(sum(f{1}.fo >= 2 & f{1}.fo <= 4)/(2./FO.resolution));
spec = {};
for I = 1:FO.nCh
    fo = [];
    s = [];
    for i = 1:f1:(size(f{I}.spec, 2) - f1)
        s = [s, mean(f{I}.spec(:, i:(i + f1 - 1), :), 2)];
        fo = [fo; mean(f{I}.fo(i:(i + f1 - 1)))];
    end
    spec{I} = s;
end

for i = 1:length(f)
    n = prctile(reshape(f{i}.spec, 1, []), [1, 99]);
    f{i}.spec(f{i}.spec < n(1)) = n(1);
    f{i}.spec(f{i}.spec > n(2)) = n(2);
end

FO.unsmoothedSpec = spec;
FO.originalFO = f{1}.fo;


for i = 1:nCh
    FO.spec{i} = log10(convWithIn(spec{i}, hanning(FO.hanningW)));
end


FO.fo = fo;
FO.to = f{1}.to;
FO.lims = [min(FO.to), max(FO.to)];

if nCh > 1
    min1 = min(min(FO.spec{1}(:, FO.fo <= FO.maxFreq)));
    max1 = max(max(FO.spec{1}(:, FO.fo <= FO.maxFreq)));
    for i = 2:nCh
        min2 = min(min(FO.spec{i}(:, FO.fo <= FO.maxFreq)));
        FO.spec{i} = FO.spec{i} - min2;
        max2 = max(max(FO.spec{i}(:, FO.fo <= FO.maxFreq)));
        FO.spec{i} = FO.spec{i}./max2;
        FO.spec{i} = FO.spec{i}.*(max1 - min1);
        FO.spec{i} = FO.spec{i} + min1;
    end
end

colors.states{1} = reshape([0, 0, 0], [1, 1, 3])/255;%black/wake
colors.states{2} = reshape([255, 236, 79], [1, 1, 3])/255;%yellow/drowzy
colors.states{3} = reshape([6, 113, 148], [1, 1, 3])/255;%blue/NREM
colors.states{4} = reshape([19, 166, 50], [1, 1, 3])/255;%green/Intermediate
colors.states{5} = reshape([207, 46, 49], [1, 1, 3])/255;%red/REM

colors.white = reshape([1, 1, 1], [1, 1, 3]);
colors.grey = reshape([206, 206, 206], [1, 1, 3])/255;
colors.orange = reshape([238, 113, 25], [1, 1, 3])/255;
FO.colors = colors;

rows{1} = 1;
rows{2} = 2;
rows{3} = 3;
rows{4} = 4;
rows{5} = 5;

FO.rows = rows;
FO.SM = ones(length(rows), length(States), 3); %SM holds the color matrix for the states plot


for I = 1:5 %sets up the states plot
    f = find(States == I);
    FO.SM(rows{I}, f, :) = repmat(colors.states{I}, [1, length(f), 1]);
end




FO.lax = axes('Position', position.lax);

FO.ilab = image(FO.to, 1:(size(FO.SM, 1)), FO.SM);
ylabel('State');
set(FO.ilab, 'HitTest', 'off');
set(FO.lax, 'YTick', [1, 2, 3, 4, 5], 'YTickLabel', (1:5), 'XTick', []);
FO.zoomL = zoom;
setAxesZoomMotion(FO.zoomL, FO.lax, 'horizontal');

FO.panL = pan;
setAxesPanMotion(FO.panL, FO.lax, 'horizontal');




for i = 1:nCh
    FO.sax{i} = axes('Position', position.sax{i});
    FO.iSpec{i} = imagesc(FO.to, FO.fo(FO.fo <= FO.maxFreq), FO.spec{i}(:, FO.fo <= FO.maxFreq)');
    colormap('jet');
    ylabel('Freq. (Hz)');
    set(FO.iSpec{i}, 'HitTest', 'off');
    set(FO.sax{i}, 'YDir', 'normal');
    FO.zoomS = zoom;
    setAxesZoomMotion(FO.zoomS, FO.sax{i}, 'horizontal');
    set(FO.zoomS, 'ActionPostCallback', {@ZoomButton});
    %set(FO.zoomS, 'Enable', 'on');
    FO.panS = pan;
    setAxesPanMotion(FO.panS, FO.sax{i}, 'horizontal');
    
    hold on;
    p = get(FO.sax{i}, 'Position');
    %    FO.sMidline{i} = annotation('line', [p(1) + p(3)/2, p(1) + p(3)/2], [p(2), p(2) + p(4)], 'LineStyle', '--', 'Color', 'w');
    FO.sMidline{i} = annotation('line', [0, 0], [p(2), p(2) + p(4)], 'LineStyle', '--', 'Color', 'w');
    set(FO.sMidline{i}, 'Position', [p(1) + p(3)/2, p(2), 0, p(4)]);
    
%     if i ~= nCh
        set(FO.sax{i}, 'XTick', [],'xticklabels',[]);
%     end
%     if i == 3
%         set(FO.sax{i}, 'XTick', [],'xticklabels',[]);
%     end
end



MP(~isnan(MP)) = zscore(MP(~isnan(MP)));
FO.max = axes('Position', position.MP);
if ~isempty(MP)
    FO.Mplot = plot(FO.to, MP, '-k');
    ylabel('Motion (z.s.)');
    set(FO.Mplot, 'HitTest', 'off');
    ylim([prctile(MP(~isnan(MP)), 1), prctile(MP(~isnan(MP)), 99)]);
    FO.zoomM = zoom;
    setAxesZoomMotion(FO.zoomM, FO.max, 'horizontal');
    FO.panM = pan;
    setAxesPanMotion(FO.panM, FO.max, 'horizontal');
    yl = get(FO.max, 'YLim');
    xl = get(FO.max, 'XLim');
    
end



p = get(FO.max, 'Position');
FO.mMidline = annotation('line', [p(1) + p(3)/2, p(1) + p(3)/2], [p(2), p(2) + p(4)], 'LineStyle', '--', 'Color', 'k');


for i = 1:nCh
    eegX = (1:length(FO.eeg{i}))/(FO.eegFS/FO.downsample);
    FO.eax{i} = axes('Position', position.eax{i});
    
    FO.Eplot{i} = plot(eegX(eegX >= 0 & eegX <= 120), FO.eeg{i}(eegX >= 0 & eegX <= 120), 'y');
    set(FO.eax{i}, 'Color', [0 0 0], 'XColor', 'b');
    %   FO.Eplot{i} = plot(eegX, FO.eeg{i});
    ylabel('Eeg');
    l1 = [min(get(FO.Eplot{i}, 'YData')), max(get(FO.Eplot{i}, 'YData'))];
    ylim(l1);
    FO.zoomE = zoom;
    setAxesZoomMotion(FO.zoomE, FO.eax{i}, 'horizontal');
    FO.panE = pan;
    setAxesPanMotion(FO.panE, FO.eax{i}, 'horizontal');
    p = get(FO.eax{i}, 'Position');
    FO.eMidline{i} = annotation('line', [p(1) + p(3)/2, p(1) + p(3)/2], [p(2), p(2) + p(4)], 'LineStyle', '--', 'Color', 'w');
    FO.lineParent = get(FO.eMidline{i}, 'Parent');
    %set(FO.lineParent, 'HandleVisibility', 'on');
    if i ~= nCh
        set(FO.eax{i}, 'XTickLabel', []);
    end
    
    if ~FO.EegUpdateBoolean
        set(FO.Eplot{i}, 'XData', [], 'YData', []);
    end
end

switch nCh
    case 1
        linkaxes([FO.lax, FO.sax{1}, FO.max, FO.lax], 'x');
    case 2
        linkaxes([FO.sax{1}, FO.sax{2}, FO.max, FO.lax], 'x');
        linkaxes([FO.eax{1}, FO.eax{2}], 'x');
    case 3
        linkaxes([FO.sax{1}, FO.sax{2}, FO.sax{3}, FO.max, FO.lax], 'x');
        linkaxes([FO.eax{1}, FO.eax{2}, FO.eax{3}], 'x');
end



%%%%% Right panel display/buttons

%Toggle for selectable states
FO.stateAxisActiveCheckbox = uicontrol('style', 'checkbox', 'String', 'Clickable State Display',...
    'value',FO.stateAxisToggleModeBool ,'Units', 'normalized', 'Position',  [0.855, 0.965, 0.1, 0.02]);
set(FO.stateAxisActiveCheckbox, 'Callback', {@ToggleStateAxisActive});

%status/action display (upper right)
FO.actionDisp = annotation('textbox', 'Units', 'normalized', 'Position', [0.85, 0.86, 0.135, 0.105], 'EdgeColor', 'none');

% Notificaton to press h for help
a = annotation('textbox', 'Position', [0.85          0.81          0.15         0.05], 'EdgeColor', 'none');
set(a, 'String', {'\fontsize{14}!Press \bf''h'' for help!';...
                    '\fontsize{14}!Press \bf''c'' to cancel'});


% Apparently unused - BW 2018
% FO.startLocDisp = annotation('textbox', 'Units', 'normalized', 'Position', [0.855, 0.7, 0.135, 0.1], 'EdgeColor', 'none');

%Guide regarding states/colors
FO.infoDisp = annotation('textbox', 'Position', [0.855, 0.7, 0.135, 0.1], 'EdgeColor', 'none');
info = {};
info{end + 1} = ['\color[rgb]{', num2str(colors.states{1}), '}1: Wake'];
info{end + 1} = ['\color[rgb]{', num2str(colors.states{2}), '}2: Drowzy/Light'];
info{end + 1} = ['\color[rgb]{', num2str(colors.states{3}), '}3: NREM'];
info{end + 1} = ['\color[rgb]{', num2str(colors.states{4}), '}4: Intermediate'];
info{end + 1} = ['\color[rgb]{', num2str(colors.states{5}), '}5: REM'];

set(FO.infoDisp, 'FontSize', 10, 'String', info);

%last click indicator
FO.lastClickDisp = annotation('textbox', 'Units', 'normalized', 'Position', [0.855, 0.68, 0.135, 0.03], 'EdgeColor', 'none');

%Go To Second Command and label
a = annotation('textbox', 'Units', 'normalized', 'Position', [0.855, 0.62, 0.135, 0.03], 'EdgeColor', 'none');
set(a, 'String', 'Go To Second:');
FO.gotosecondbox = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.88, 0.605, 0.06, 0.025]);
set(FO.gotosecondbox, 'Callback', @goToSecond);

% Window length command and labe
a = annotation('textbox', 'Units', 'normalized', 'Position', [0.855, 0.57, 0.135, 0.03], 'EdgeColor', 'none');
set(a, 'String', 'Window Length (sec)');
FO.xlimbox = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.88, 0.55, 0.06, 0.025]);
set(FO.xlimbox, 'Callback', {@changeXlim}, 'String', int2str(round(diff(FO.lims))));

%Smoothing windows settings stuff
Woptions =  [0, 5, 10, 15, 20, 30, 45, 60];
FO.Woptions = Woptions;
optString = [];
for I = 1:length(Woptions)
    optString = [optString, int2str(Woptions(I)),' secs|'];
end
optString = optString(1:(end - 1));
hanL = annotation('textbox', 'Units', 'normalized', 'Position', [0.855, 0.515, 0.135, 0.03], 'EdgeColor', 'none','String', 'Smoothing Window:');
FO.hanningWDisp = uicontrol('style', 'popup', 'Units', 'normalized', 'Position', [0.88, 0.505, 0.08, 0.01]);
set(FO.hanningWDisp, 'String', optString, 'CallBack', {@ChangeSmoothingWindow}, 'Value', find(Woptions == FO.hanningW));

%Overlay stuff
Ooptions = ['none|(5-10Hz)/(0.5-4Hz)|From SleepScoreMaster|Choose from file'];
a = annotation('textbox', 'Units', 'normalized', 'Position', [0.855, 0.46, 0.1355, 0.03], 'EdgeColor', 'none');
set(a, 'String', 'Overlay Display:');
FO.overlayDisp = uicontrol('style', 'popup', 'Units', 'normalized', 'Position', [0.8800    0.45    0.0800    0.01]);
set(FO.overlayDisp, 'String', Ooptions, 'CallBack', {@OverlayDisplay}, 'Value', 1);

%Event stuff
Eoptions = ['none|1 (0 events)|2 (0 events)|3 (0 events)|4 (0 events)|5 (0 events)|6 (0 events)|7 (0 events)|8 (0 events)|9 (0 events)|10 (0 events)'];
a = annotation('textbox', 'Units', 'normalized', 'Position', [0.855, 0.403, 0.1355, 0.03], 'EdgeColor', 'none');
set(a, 'String', 'Event #:');
FO.eventDisp = uicontrol('style', 'popup', 'Units', 'normalized', 'Position', [0.8800    0.397    0.0800    0.01]);
set(FO.eventDisp, 'String', Eoptions, 'CallBack', {@EventNumber}, 'Value', 2);
FO.eventNum = 1;

% EEG stuff
%Set width of eeg display 
a = annotation('textbox', 'Units', 'normalized', 'Position', [0.87, 0.225, 0.135, 0.03], 'EdgeColor', 'none');
set(a, 'String', 'EEG Display Width (sec):');
FO.eegDisplayWidthBox = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.88, 0.21, 0.06, 0.025]);
set(FO.eegDisplayWidthBox,'String',num2str(FO.eegDisplaySeconds))
set(FO.eegDisplayWidthBox, 'Callback', {@changeEEGDisplaySeconds});

%Toggle for whether or not to update eeg when clicking
FO.updateEegOnClickCheckbox = uicontrol('style', 'checkbox', 'String', 'Update EEG on Click', 'Units', 'normalized', 'Position',  [0.87, 0.185, 0.1, 0.02]);
set(FO.updateEegOnClickCheckbox, 'Callback', {@ToggleUpdateEeg});

%Save raw eeg for offline use
FO.saveEEGButton = uicontrol('style', 'pushbutton', 'String', 'Export Raw EEG', 'Units', 'normalized', 'Position',  [0.87, 0.152, 0.1, 0.03]);
set(FO.saveEEGButton, 'Callback', {@saveRawEEG});



% Buttons to scroll full screen widths
a = annotation('textbox', 'Units', 'normalized', 'Position', [0.87, 0.078, 0.2, 0.03], 'EdgeColor', 'none');
set(a, 'String', 'Scroll:');
FO.leftScreenOver = uicontrol('style', 'pushbutton', 'String', '<', 'Units', 'normalized', 'Position',  [0.898, 0.08, 0.038, 0.03]);
set(FO.leftScreenOver, 'Callback', {@leftScreenOver});
FO.rightScreenOver = uicontrol('style', 'pushbutton', 'String', '>', 'Units', 'normalized', 'Position',  [0.938, 0.08, 0.038, 0.03]);
set(FO.rightScreenOver, 'Callback', {@rightScreenOver});

%Undo/Redo state buttons
FO.undoButton = uicontrol('style', 'pushbutton', 'String', 'Undo State', 'Units', 'normalized', 'Position',  [0.87, 0.05, 0.1, 0.03]);
set(FO.undoButton, 'Callback', {@undoChange});
FO.redoButton = uicontrol('style', 'pushbutton', 'String', 'Redo State', 'Units', 'normalized', 'Position',  [0.87, 0.02, 0.1, 0.03]);
set(FO.redoButton, 'Callback', {@redoChange});

%???
for i = 1:length(FO.Chs)
    a = annotation('textbox', 'Units', 'Normalized');
    if iscell(FO.Chs)
        set(a, 'String', ['\bf\color{red}\fontsize{10}C', FO.Chs{i}], 'Position', position.eegCh{i}, 'EdgeColor', 'none');        
    else
        set(a, 'String', ['\bf\color{red}\fontsize{10}C', int2str(FO.Chs(i))], 'Position', position.eegCh{i}, 'EdgeColor', 'none');
    end
end
FO.eegWidthDisp = annotation('textbox', 'Units', 'Normalized');
set(FO.eegWidthDisp, 'String', ['\bf\color{red}\fontsize{11}', num2str(FO.eegDisplaySeconds), ' sec'], 'Position', position.eegWidth, 'EdgeColor', 'none');

%%%% Other stuff?

%set(FO.actionDisplay, 'String', {'\fontsize{12}\bfCurrent Action:', ' ', '\fontsize{20}Browse'});

FO.Events = [];
FO.CurrEventLines = {};
FO.Transitions = [];
FO.TransHistoryTracker = [];
FO.saxYLim = [min(FO.fo), max(FO.fo)];
FO.mpYLim = get(FO.max, 'YLim');
a = [];
for i = 1:length(nCh)
    a = [a; min(FO.eeg{i}), max(FO.eeg{i})];
end
FO.eegYLim = [min(a(:, 1)), max(a(:, 2))];


%% BW speeding things up... didn't change anything above to be safe
setappdata(FO.fig,'unsmoothedSpec',FO.unsmoothedSpec)
setappdata(FO.fig,'spec',FO.spec)
setappdata(FO.fig,'eeg',FO.eeg)

FO = rmfield(FO,'spec');
FO = rmfield(FO,'unsmoothedSpec');%access only when needed using appdata now
FO = rmfield(FO,'eeg');%access only when needed using appdata now
% FO = rmfield(FO,'eegX');%recalculate on the fly using:   eegX = (1:length(FO.eeg{i}))/(FO.eegFS/FO.downsample);
%%
guidata(FO.fig, FO); 


if FO.EegUpdateBoolean
    FO = updateEEG(FO);
end
FO = UpdateGUI(FO);
set(FO.max,'xticklabel',num2str(get(FO.max,'xtick')'));
% set(FO.sax{end},'xticklabel',num2str(get(FO.sax{end},'xtick')'));
set(FO.sax{1}, 'XLim', FO.lims);

%set(FO.eax{end},'xticklabel',num2str(get(FO.eax{end},'xtick')'));



end


function DefKey(fig, e)
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
FO = guidata(fig);

switch e.Key
    case 'uparrow'
        if strcmpi(FO.currAction, 'FreqResize')
            ResizeFreqY(1);
            obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
        else
            for i = 1:FO.nCh
                v = caxis(gca); caxis(gca, v - 0.1);
            end
        end
    case 'downarrow'
        if strcmpi(FO.currAction, 'FreqResize')
            ResizeFreqY(-1);
            obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
        else
            for i = 1:FO.nCh
                v = caxis(gca); caxis(gca, v + 0.1);
            end
        end
    case 'z'
        if strcmpi(FO.currAction, 'Browse');
            for i = 1:FO.nCh
                axes(FO.sax{i});
            end
            %zoom;
            FO.currAction = 'Zoom';
            
        else
            for i = 1:FO.nCh
                axes(FO.sax{i});
                
            end
            
            if isempty(FO.startLocation)
                FO.currAction = 'Browse';
            else
                FO.currAction = 'Add';
            end
        end
    case 'rightarrow'
        l = get(FO.sax{1}, 'XLim');
        l = l + 0.15*diff(l);
        if l(1) < FO.lims(1)
            l(1) = FO.lims(1);
        end
        if l(2) > FO.lims(2)
            l(2) = FO.lims(2);
        end
        set(FO.sax{1}, 'XLim', l);
%         guidata(FO.fig, FO); 
        if FO.EegUpdateBoolean
            FO = updateEEG(FO);
        end
%         obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
        % updateMidline([]);
    case 'leftarrow'
        l = get(FO.sax{1}, 'XLim');
        l = l - 0.15*diff(l);
        if l(1) < FO.lims(1)
            l(1) = FO.lims(1);
        end
        if l(2) > FO.lims(2)
            l(2) = FO.lims(2);
        end
        set(FO.sax{1}, 'XLim', l);
%         guidata(FO.fig, FO); 
        if FO.EegUpdateBoolean
            FO = updateEEG(FO,nextE);
        end        
%         obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
        %    updateMidline([]);
    case '1'
        FO.currentState = 1;
        FO.currAction = 'Add';
        FO.stateAxisToggleModeBool = logical(0);
        
        set(FO.fig,'Pointer','arrow');
        FO.zoomState = 0;
        set(FO.fig, 'Name', ['States: ', FO.baseName, ' - Add State 1']);
    case '2'
        FO.currentState = 2;
        FO.currAction = 'Add';
        FO.stateAxisToggleModeBool = logical(0);
        
        set(FO.fig,'Pointer','arrow');
        set(FO.fig, 'Name', ['States: ', FO.baseName, ' - Add State 2']);
    case '3'
        FO.currentState = 3;
        FO.currAction = 'Add';
        FO.stateAxisToggleModeBool = logical(0);
        
        set(FO.fig,'Pointer','arrow');
        FO.zoomState = 0;
        set(FO.fig, 'Name', ['States: ', FO.baseName, ' - Add State 3']);
    case '4'
        FO.currentState = 4;
        FO.currAction = 'Add';
        FO.stateAxisToggleModeBool = logical(0);
        
        set(FO.fig,'Pointer','arrow');
        FO.zoomState = 0;
        set(FO.fig, 'Name', ['States: ', FO.baseName, ' - Add State 4']);
    case '5'
        FO.currentState = 5;
        FO.currAction = 'Add';
        FO.stateAxisToggleModeBool = logical(0);
        
        set(FO.fig,'Pointer','arrow');
        FO.zoomState = 0;
        set(FO.fig, 'Name', ['States: ', FO.baseName, ' - Add State 5']);
    case '0'
        FO.currentState = 0;
        FO.currAction = 'Add';
        
        set(FO.fig,'Pointer','arrow');
        FO.zoomState = 0;
        set(FO.fig, 'Name', ['States: ', FO.baseName, ' - Add State 0']);
    case 'e'
        FO.currentState = 0;
        if strcmp(FO.currAction, 'AddEvent')
            if isempty(FO.startLocation)
                FO.currAction = 'Browse';
            else
                FO.currAction = 'Add';
            end
            
        else
            if length(FO.eventNum) == 1
                FO.currAction = 'AddEvent';
                set(gcf, 'Pointer', 'crosshair');
            else
                warndlg('No event # selected');
            end
        end
        FO = ResetStateAxisClicabilityByCheckboxStatus(FO);
    case 'd'
        if strcmp(FO.currAction, 'DeleteEvent')
            if isempty(FO.startLocation)
                FO.currAction = 'Browse';
            else
                FO.currAction = 'Add';
            end
        else
            if ~isempty(FO.Events)
                if sum(FO.Events(:, 1) == FO.eventNum) > 0
                    skullCursor;
                    FO.currAction = 'DeleteEvent';
                else
                    warndlg('There are currently no events displayed for you to delete');
                end
            else
                warndlg('There are currently no events displayed for you to delete');
            end
        end
        FO = ResetStateAxisClicabilityByCheckboxStatus(FO);

    case {'c', 'escape'}
        FO.currAction = 'Browse';
        if ~isempty(FO.startLine)
            for i = 1:length(FO.startLine)
                delete(FO.startLine{i});
            end
            FO.startLine = {};
            
        end
        FO = ResetStateAxisClicabilityByCheckboxStatus(FO);
    case 'a'
        FO = ViewAutoScoreThresholds(gcf);
    case 'f'
        if strcmpi(FO.currAction, 'FreqResize')
            FO.currAction = 'Browse';
        else
            FO.currAction = 'FreqResize';
        end
    case 'u'
        undoChange;
        obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
        FO = ResetStateAxisClicabilityByCheckboxStatus(FO);
    case 'r'
        set(FO.sax{1}, 'XLim', FO.lims);
    case 's'
        saveStates;
        obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
    case 'l'
        LoadStates;
        obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
        FO = ResetStateAxisClicabilityByCheckboxStatus(FO);
    case 'n'
        nextEvent;
    case 'p'
        previousEvent;
    case 'hyphen'
        f = find(histc(FO.eegDisplaySeconds, [0, 0.26, 2.1, 5.1, 15.1, 30.1,  60.1]) == 1);
        switch f
            case 1
                delta = 0;
            case 2
                delta = 0.25;
            case 3
                delta = 0.5;
            case 4
                delta = 1;
            case 5
                delta = 2.5;
            case 6
                delta = 5;
        end
        
        
        FO.eegDisplaySeconds = FO.eegDisplaySeconds - delta;
        guidata(FO.fig, FO); 
        if FO.EegUpdateBoolean
            FO = updateEEG(FO,mean(get(FO.eax{1}, 'XLim')));
        end
        FO = UpdateGUI(FO);
        return;
        
    case 'equal'
        f = find(histc(FO.eegDisplaySeconds, [0, 0.24, 1.99, 4.99, 14.99, 29.9,  60]) == 1);
        switch f
            case 2
                delta = 0.25;
            case 3
                delta = 0.5;
            case 4
                delta = 1;
            case 5
                delta = 2.5;
            case 6
                delta = 5;
            case 7
                delta = 0;
        end 
        
        FO.eegDisplaySeconds = FO.eegDisplaySeconds + delta;
        guidata(FO.fig, FO); 
        
        if FO.EegUpdateBoolean
            FO = updateEEG(FO,mean(get(FO.eax{1}, 'XLim')));
        end
        FO = UpdateGUI(FO);;
        return;
        
    case 'h'
        helpFig = figure('Position', [300    50   520   750], 'MenuBar', 'none');
        a = annotation('textbox', 'Units', 'normalized');
        set(a, 'Position', [0.025 0.05 0.95 0.95], 'EdgeColor', 'none');
        set(a, 'String', {'\fontsize{14}\bfWelcome to StateEditor!\fontsize{10}\rm'...
            ' '...
            '''\bfSingle-Click\rm''-- Center view on click point'...
            ' ',...
            '''\bfClick-and-Hold\rm''-- Drag currently visible ''xlim'' extent(applied to either spectral'...
            '        or lfp windows depending on mouse position).'...
            ' '...
            '''\bfScroll-Wheel\rm''-- Zoom in or out (applied to either spectral or lfp windows depending'...
            '        on click position).'...
            ' '...
            '''\bfLeftArrow\rm''-- Move view left'...
            '''\bfRightArrow\rm''-- Move view right'...
            ' '...
            '''\bfUpArrow\rm''-- Increase color limits'...
            '''\bfDownArrow\rm''-- Decrease color limits'...
            ' '...
            '\bf''1'', ''2'', ''3'', ''4'' or ''5''\rm',...
            '       Add State - First click adds first bound, second click adds second bound.'...
            '       Press ''\bfC'' \rmto cancel addition.'...
            ' '...
            '\bf''0''\rm--Delete state labels (add state 0)'...
            ' '...
            '''\bfZ''\rm-- Toggle Zoom ON/OFF',...
            '       Left click zoom in. Right click zoom out. Hold and drag to select zoom area.'...
            '       Double left click: fast zoom in. Double right click: reset full X extent.'... 
            ' '...
            '''\bfR''\rm-- Reset X limits to full extent',...
            ' '...
            '''\bfF\rm''-- Toggle frequency edit mode. ''UpArrow''- Increase Freq. extent'...
            '       ''DownArrow''- Decrease Freq. extent'...
            ' '...
            '''\bf- or =\rm''-- Increase/decrease the extent of LFP display(s).'...
            ' '...
            '''\bfE\rm''-- Add event (the number is determined by the Event # list on the right'...
            '       panel.'...
            ' '...
            '''\bfD\rm''-- Delete event from currently selected event number.'...
            ' '...
            '''\bfN or P\rm''-- View next/Previous event (of currently selected event #)'...
            ' '...
            '''\bfS''\rm-- Save state vector/events/transitions to file',...
            ' '...
            '''\bfL''\rm-- Load state vector/events/transitions from file',...
            ' '...
            '''\bfA''\rm-- Autoscoring parameters: View and manipulate',...
            ' '...
            'Also note the editable fields and lists on the right hand panel!'
                
            });
        uiwait(helpFig);
        obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
end
guidata(FO.fig, FO);
FO = UpdateGUI(FO);
end



function FO = updateEEG(FO,varargin)
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
% try
    eeg = getappdata(FO.fig,'eeg');
    eegX = (1:length(eeg{1}))/(FO.eegFS/FO.downsample);
% catch
%     eeg = FO.eeg;
%     eegX = FO
% end

if isempty(varargin)
    pos = mean(get(FO.sax{1}, 'XLim'));
else
    pos = varargin{1};
end
low = pos - FO.eegDisplaySeconds/2;
high = pos + FO.eegDisplaySeconds/2;
if low < FO.lims(1)
    high = high + (FO.lims(1) + low);
    low = FO.lims(1);
else
    if high > FO.lims(2)
        low = low + (FO.lims(2) - high);
        high = FO.lims(2);
    end
end

lowMargin = low - 60;
highMargin = high + 60;
for i = 1:FO.nCh
    set(FO.Eplot{i}, 'XData', eegX(eegX >= lowMargin & eegX <= highMargin), 'YData', eeg{i}(eegX >= lowMargin & eegX <= highMargin));
    l1 = [min(eeg{i}(eegX >= low & eegX <= high)), max(eeg{i}(eegX >= low & eegX <= high))];
    set(FO.eax{i}, 'YLim', l1);
    set(FO.eax{i}, 'XLim', [pos - FO.eegDisplaySeconds/2, pos + FO.eegDisplaySeconds/2]);
end

% set(FO.max,'xticklabel',num2str(get(FO.max,'xtick')'));
% set(FO.sax{end},'xticklabel',num2str(get(FO.sax{end},'xtick')'));
set(FO.eax{end},'xticklabel',num2str(get(FO.eax{end},'xtick')'));

lims1 = get(FO.sax{1}, 'XLim');
perc = (pos - lims1(1))/diff(lims1);
newL = (diff(FO.xplotLims)*perc) + FO.xplotLims(1);

for i = 1:length(FO.sMidline)
    set(FO.sMidline{i}, 'X', [newL, newL]);
end
set(FO.mMidline, 'X', [newL, newL]);


h = get(gcf, 'Children');
try
    h = [FO.lineParent; h(h ~= FO.lineParent)];
    warning('off', 'MATLAB:hg:default_child_strategy:IllegalPermutation')
    set(gcf, 'Children', h);
catch
end
set(FO.eegDisplayWidthBox,'String',FO.eegDisplaySeconds);%update display



% for i = 1:FO.nCh
%     set(FO.Eplot{i}, 'XData', eegX(eegX >= low & eegX <= high));
%     set(FO.Eplot{i}, 'YData', eeg{i}(eegX >= low & eegX <= high));
%     l1 = [min(eeg{i}(eegX >= low & eegX <= high)), max(eeg{i}(eegX >= low & eegX <= high))];
%     set(FO.eax{i}, 'YLim', l1);
%     set(FO.eax{i}, 'XLim', [pos - FO.eegDisplaySeconds/2, pos + FO.eegDisplaySeconds/2]);
% end

% guidata(FO.fig, FO); 
end


function MouseClick(obj, ev)
st = dbstack;%make sure not overloading thestateeditor
if length(st)>1
    for idx = 2:length(st)
        if strmatch (st(idx).file,'TheStateEditor.m')%if state editor is already running
            disp('You clicked while TheStateEditor was still running, click later')
            return
        end
    end
end

% obj = findobj('tag','StateEditorMaster');  
FO = guidata(obj); 
global isClicking
persistent chk lastButton

if ~strcmpi(get(FO.fig, 'SelectionType'), 'open')
    lastButton = get(FO.fig, 'SelectionType');
end
% updateEegToClick = 1;
isClicking = 1;

clickType = [];
holdC = get(FO.fig, 'CurrentPoint');
if isempty(chk)
    chk = 1;
    pause(0.25);
    if chk == 1
        if isClicking == 0
            
            clickType = 'Single';
            
        else
            clickType = 'Hold';
        end
        
        chk = [];
    end
else
    chk = [];
    clickType = 'Double';
end
if strcmpi(clickType, 'Double')
    return;
end

if isempty(clickType)
    clickType = 'Double';
end


c = get(FO.fig, 'CurrentPoint');

clickloc = [];
if FO.stateAxisToggleModeBool
    xy = FO.lax.Position;
    xs = [xy(1) xy(1) xy(1)+xy(3) xy(1)+xy(3)];
    ys = [xy(2) xy(2)+xy(4) xy(2)+xy(4) xy(2)];
    inp = inpolygon(c(1),c(2),xs,ys);%if in top axes
    if inp
        clickloc = 'activestateaxes';
    end
end
    %get coordinate in that axes
if isempty(clickloc)
    xG = (c(1) >= FO.xplotLims(1)) & (c(1) <=  FO.xplotLims(2));
    yG = sum((c(2) >= FO.yplotLims(:, 1)) & c(2) <= FO.yplotLims(:, 2)) == 1;
    clickloc = 'outside';
    if ~(xG && yG)

        yG = sum((c(2) >= FO.yplotLFPLims(:, 1)) & c(2) <= FO.yplotLFPLims(:, 2)) ~= 0;
        if ~(xG && yG)
            return;
        else
            clickloc = 'timeaxes';
        end
    end
end

% %taking care of weird error I can't figure out -BW
% if strcmp('browse',lower(FO.currAction))
% %     if ~isempty(strfind(get(FO.fig,'name'),'Add State'))
%     if ~isempty(FO.startLine)
%         FO.startLine = {};
%     end
% end


switch(clickloc)
    case 'outside'
        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
        xl = get(FO.sax{1}, 'XLim');
        pointTo = xl(1) + (diff(xl)*d);
        button = lastButton;
        
        
        switch(FO.currAction);
            case 'Add'
                switch(clickType)
                    case 'Single'
                        FO = addStateLine(FO,pointTo);
                    case 'Double'
                        FO = addStateLine(FO,pointTo);
                    case 'Hold'
                        c = holdC;
                        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                        xl = get(FO.sax{1}, 'XLim');
                        pointTo = xl(1) + (diff(xl)*d);
                        origin = xl;
                        while isClicking == 1
                            set(gcf, 'Pointer', 'fleur');
                            c = get(gcf, 'CurrentPoint');
                            d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                            newPoint = xl(1) + (diff(xl)*d);
                            delta = pointTo - newPoint;
                            
                            if ((origin(1) + delta) >= FO.lims(1)) & ((origin(2) + delta) <= FO.lims(2))
                                set(FO.sax{1}, 'XLim', [origin(1) + delta, origin(2) + delta]);
                            end
                            set(FO.max,'xticklabel',num2str(get(FO.max,'xtick')'));
%                             set(FO.sax{end},'xticklabel',num2str(get(FO.sax{end},'xtick')'));
                            pause(0.025);
                        end
                        set(gcf, 'Pointer', 'arrow');
                        if FO.EegUpdateBoolean
                            FO = updateEEG(FO,mean(get(FO.sax{1}, 'XLim')));
                        end
%                         updateEegToClick = 0;
                end
            case 'AddEvent'
                switch(clickType)
                    case 'Single'
                        FO = addEvent(FO,pointTo);
                    case 'Double'
                        FO = addEvent(FO,pointTo);
                    case 'Hold'
                        c = holdC;
                        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                        xl = get(FO.sax{1}, 'XLim');
                        pointTo = xl(1) + (diff(xl)*d);
                        origin = xl;
                        while isClicking == 1
                            set(gcf, 'Pointer', 'fleur');
                            c = get(gcf, 'CurrentPoint');
                            d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                            newPoint = xl(1) + (diff(xl)*d);
                            delta = pointTo - newPoint;
                            
                            if ((origin(1) + delta) >= FO.lims(1)) & ((origin(2) + delta) <= FO.lims(2))
                                set(FO.sax{1}, 'XLim', [origin(1) + delta, origin(2) + delta]);
                            end
                            set(FO.max,'xticklabel',num2str(get(FO.max,'xtick')'));
%                             set(FO.sax{end},'xticklabel',num2str(get(FO.sax{end},'xtick')'));
                            pause(0.025);
                        end
                        set(gcf, 'Pointer', 'crosshair');
                        if FO.EegUpdateBoolean
                            FO = updateEEG(FO,mean(get(FO.sax{1}, 'XLim')));
                        end
%                         updateEegToClick = 0;
                end
            case 'DeleteEvent'
                switch(clickType)
                    case 'Single'
                        FO = deleteEvent(FO,pointTo);
                    case 'Double'
                        FO = deleteEvent(FO,pointTo);
                    case 'Hold'
                        c = holdC;
                        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                        xl = get(FO.sax{1}, 'XLim');
                        pointTo = xl(1) + (diff(xl)*d);
                        origin = xl;
                        while isClicking == 1
                            set(gcf, 'Pointer', 'fleur');
                            c = get(gcf, 'CurrentPoint');
                            d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                            newPoint = xl(1) + (diff(xl)*d);
                            delta = pointTo - newPoint;
                            
                            if ((origin(1) + delta) >= FO.lims(1)) & ((origin(2) + delta) <= FO.lims(2))
                                set(FO.sax{1}, 'XLim', [origin(1) + delta, origin(2) + delta]);
                            end
                            set(FO.max,'xticklabel',num2str(get(FO.max,'xtick')'));
%                             set(FO.sax{end},'xticklabel',num2str(get(FO.sax{end},'xtick')'));
                            pause(0.025);
                        end
                        set(gcf, 'Pointer', 'crosshair');
                        if FO.EegUpdateBoolean
                            FO = updateEEG(FO,mean(get(FO.sax{1}, 'XLim')));
                        end
%                         updateEegToClick = 0;
                end
            case 'Zoom'
                switch(clickType)
                    case 'Single'
                        switch(button)
                            case 'alt'
                                xl = xl*1.25;
                                X1 = ((pointTo) - diff(xl)/2);
                                X2 = ((pointTo) + diff(xl)/2);
                                if X1 < FO.lims(1)
                                    X2 = X2 + (FO.lims(1) - X1);
                                    X1 = FO.lims(1);
                                end
                                
                                if X2 > FO.lims(2)
                                    X1 = X1 - (X2 - FO.lims(2));
                                    X2 = FO.lims(2);
                                end
                                if X1 < FO.lims(1)
                                    X1 = FO.lims(1);
                                end
                                set(FO.sax{1}, 'XLim', [X1, X2]);
                            case 'normal'
                                xl = xl*0.75;
                                X1 = ((pointTo) - diff(xl)/2);
                                X2 = ((pointTo) + diff(xl)/2);
                                if X1 < FO.lims(1)
                                    X2 = X2 + (FO.lims(1) - X1);
                                    X1 = FO.lims(1);
                                end
                                
                                if X2 > FO.lims(2)
                                    X1 = X1 - (X2 - FO.lims(2));
                                    X2 = FO.lims(2);
                                end
                                
                                if X1 < FO.lims(1)
                                    X1 = FO.lims(1);
                                end
                                set(FO.sax{1}, 'XLim', [X1, X2]);
                        end
                    case 'Double'
                        switch(button)
                            case 'alt'
                                set(FO.sax{1}, 'XLim', FO.lims);
                            case 'normal'
                                xl = xl*0.25;
                                X1 = ((pointTo) - diff(xl)/2);
                                X2 = ((pointTo) + diff(xl)/2);
                                if X1 < FO.lims(1)
                                    X2 = X2 + (FO.lims(1) - X1);
                                    X1 = FO.lims(1);
                                end
                                
                                if X2 > FO.lims(2)
                                    X1 = X1 - (X2 - FO.lims(2));
                                    X2 = FO.lims(2);
                                end
                                
                                if X1 < FO.lims(1)
                                    X1 = FO.lims(1);
                                end
                                set(FO.sax{1}, 'XLim', [X1, X2]);
                        end
                    case 'Hold'
                        switch(button)
                            case 'normal'
                                if FO.EegUpdateBoolean
                                    FO = updateEEG(FO);
                                end
                                set(gcf, 'Pointer', 'circle');
                                point1 = get(gcf,'CurrentPoint');    % button down detected
%                                 finalRect = rbbox;                   % return figure units
                                point2 = get(gcf,'CurrentPoint');    % button up detected
                                point1 = point1(1,1);
                                point2 = point2(1,1);
                                d1 = (point1 - FO.xplotLims(1))./diff(FO.xplotLims);
                                d2 = (point2 - FO.xplotLims(1))./diff(FO.xplotLims);
                                xl = get(FO.sax{1}, 'XLim');
                                pointTo1 = xl(1) + (diff(xl)*d1);
                                pointTo2 = xl(1) + (diff(xl)*d2);
                                if pointTo1 > pointTo2
                                    p = pointTo2;
                                    pointTo2 = pointTo1;
                                    pointTo1 = p;
                                end
                                if pointTo1 < FO.lims(1)
                                    pointTo1 = FO.lims(1);
                                end
                                if pointTo2 > FO.lims(2)
                                    pointTo2 = FO.lims(2);
                                end
                                set(gcf, 'Pointer', 'crosshair');
                                set(FO.sax{1}, 'XLim', [pointTo1, pointTo2]);
                        end
                end
                
            case 'Browse'
                if strcmp(clickType, 'Double')
                    clickType = 'Single';
                end
                switch(clickType)
                    case 'Single'
                        if (((pointTo) - diff(xl)/2) >= FO.lims(1)) & (((pointTo) + diff(xl)/2) <= FO.lims(2))
                            set(FO.sax{1}, 'XLim', [((pointTo) - diff(xl)/2), ((pointTo) + diff(xl)/2)]);
                        end
                    case 'Hold'
                        if FO.EegUpdateBoolean
                            FO = updateEEG(FO);
                        end
                        c = holdC;
                        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                        xl = get(FO.sax{1}, 'XLim');
                        pointTo = xl(1) + (diff(xl)*d);
                        origin = xl;
                        while isClicking == 1
                            set(gcf, 'Pointer', 'fleur');
                            c = get(gcf, 'CurrentPoint');
                            d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                            newPoint = xl(1) + (diff(xl)*d);
                            delta = pointTo - newPoint;
                            set(FO.max,'xticklabel',num2str(get(FO.max,'xtick')'));
%                             set(FO.sax{end},'xticklabel',num2str(get(FO.sax{end},'xtick')'));
                            if ((origin(1) + delta) >= FO.lims(1)) & ((origin(2) + delta) <= FO.lims(2))
                                set(FO.sax{1}, 'XLim', [origin(1) + delta, origin(2) + delta]);
                            end
                            set(FO.max,'xticklabel',num2str(get(FO.max,'xtick')'));
%                             set(FO.sax{end},'xticklabel',num2str(get(FO.sax{end},'xtick')'));
                            
                            pause(0.01);
                        end
                        set(gcf, 'Pointer', 'hand');
                        if FO.EegUpdateBoolean
                            FO = updateEEG(FO,mean(get(FO.sax{1}, 'XLim')));
                        end
%                         updateEegToClick = 0;
                end
        end
    case 'timeaxes'
        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
        xl = get(FO.eax{1}, 'XLim');
        pointTo = xl(1) + (diff(xl)*d);
        button = lastButton;
        switch(FO.currAction);
            case 'Add'
                switch(clickType)
                    case 'Single'
                        FO = addStateLine(FO,pointTo);
                    case 'Hold'
                        c = holdC;
                        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                        xl = get(FO.eax{1}, 'XLim');
                        pointTo = xl(1) + (diff(xl)*d);
                        origin = xl;
                        while isClicking == 1
                            set(gcf, 'Pointer', 'fleur');
                            c = get(gcf, 'CurrentPoint');
                            d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                            newPoint = xl(1) + (diff(xl)*d);
                            delta = pointTo - newPoint;
                            
                            if ((origin(1) + delta) >= FO.lims(1)) & ((origin(2) + delta) <= FO.lims(2))
                                set(FO.eax{1}, 'XLim', [origin(1) + delta, origin(2) + delta]);
                            end
                            pause(0.025);
                        end
                        set(gcf, 'Pointer', 'arrow');
                end
            case 'AddEvent'
                switch(clickType)
                    case 'Single'
                        FO = addEvent(FO,pointTo);
                    case 'Double'
                        FO = addEvent(FO,pointTo);
                    case 'Hold'
                        c = holdC;
                        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                        xl = get(FO.eax{1}, 'XLim');
                        pointTo = xl(1) + (diff(xl)*d);
                        origin = xl;
                        while isClicking == 1
                            set(gcf, 'Pointer', 'fleur');
                            c = get(gcf, 'CurrentPoint');
                            d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                            newPoint = xl(1) + (diff(xl)*d);
                            %    updateEEG(newPoint);
                            delta = pointTo - newPoint;
                            
                            if ((origin(1) + delta) >= FO.lims(1)) & ((origin(2) + delta) <= FO.lims(2))
                                set(FO.eax{1}, 'XLim', [origin(1) + delta, origin(2) + delta]);
                            end
                            set(FO.eax{end},'xticklabel',num2str(get(FO.eax{end},'xtick')'));
                            pause(0.025);
                        end
                        
                        set(gcf, 'Pointer', 'crosshair');
                end
            case 'DeleteEvent'
                switch(clickType)
                    case 'Single'
                        FO = deleteEvent(pointTo);
                    case 'Double'
                        FO = deleteEvent(pointTo);
                    case 'Hold'
                        c = holdC;
                        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                        xl = get(FO.eax{1}, 'XLim');
                        pointTo = xl(1) + (diff(xl)*d);
                        origin = xl;
                        while isClicking == 1
                            set(gcf, 'Pointer', 'fleur');
                            c = get(gcf, 'CurrentPoint');
                            d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                            newPoint = xl(1) + (diff(xl)*d);
                            %    updateEEG(newPoint);
                            delta = pointTo - newPoint;
                            
                            if ((origin(1) + delta) >= FO.lims(1)) & ((origin(2) + delta) <= FO.lims(2))
                                set(FO.eax{1}, 'XLim', [origin(1) + delta, origin(2) + delta]);
                            end
                            set(FO.eax{end},'xticklabel',num2str(get(FO.eax{end},'xtick')'));
                            pause(0.025);
                        end
                        
                        skullCursor;
                end
            case 'Browse'
                if strcmp(clickType, 'Double')
                    clickType = 'Single';
                end
                switch(clickType)
                    case 'Single'
                        if (((pointTo) - diff(xl)/2) >= FO.lims(1)) & (((pointTo) + diff(xl)/2) <= FO.lims(2))
                            set(FO.eax{1}, 'XLim', [((pointTo) - diff(xl)/2), ((pointTo) + diff(xl)/2)]);
                        end
                    case 'Hold'
                        c = holdC;
                        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                        xl = get(FO.eax{1}, 'XLim');
                        pointTo = xl(1) + (diff(xl)*d);
                        origin = xl;
                        while isClicking == 1
                            set(gcf, 'Pointer', 'fleur');
                            c = get(gcf, 'CurrentPoint');
                            d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
                            newPoint = xl(1) + (diff(xl)*d);
                            %    updateEEG(newPoint);
                            delta = pointTo - newPoint;
                            
                            if ((origin(1) + delta) >= FO.lims(1)) & ((origin(2) + delta) <= FO.lims(2))
                                set(FO.eax{1}, 'XLim', [origin(1) + delta, origin(2) + delta]);
                            end
                            set(FO.eax{end},'xticklabel',num2str(get(FO.eax{end},'xtick')'));
                            pause(0.025);
                        end
                        
                        set(gcf, 'Pointer', 'hand');
                end
        end
        pointTo = mean(get(FO.eax{1}, 'XLim'));
        if FO.EegUpdateBoolean
            FO = updateEEG(FO,pointTo);
        end
        xl = (get(FO.sax{1}, 'XLim'));
        if (((pointTo) - diff(xl)/2) >= FO.lims(1)) & (((pointTo) + diff(xl)/2) <= FO.lims(2))
            set(FO.sax{1}, 'XLim', [((pointTo) - diff(xl)/2), ((pointTo) + diff(xl)/2)]);
        end
    case 'activestateaxes'
        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
        xl = get(FO.lax, 'XLim');
        pointTo = xl(1) + (diff(xl)*d);
        pointTo = round(pointTo);
        
        FO = ChangeThisStateAssigmnent_Setup(FO,pointTo);
%         states = FO.States;
%         ds = logical(abs(diff(states)));%find state switches with 1's
%         ds = cat(2,ds,0);%pad
%         %find start
%         spanstart = find(ds(1:pointTo)>0,1,'last')+1;
%         if isempty(spanstart)
%             spanstart = 1;
%         end
%         if spanstart<1
%             spanstart = 1;
%         end
%         %find end
%         spanend = find(ds(pointTo:end)>0,1,'first') + pointTo-1;
%         if isempty(spanend)
%             spanend = length(states);;
%         end
%         if spanend>length(states)
%             spanend = length(states);
%         end
% 
%         set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', ' ', '\fontsize{20}Choose New State (1-5)'});
%         
%         set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', ' ', '\fontsize{20}Browse'});

        
%         button = lastButton;
        % execute state switching:
            %find span occupying this click
            %ask user what to change it to with 1-5 input
            %if anything else just run the defkey='c' code

end
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
% disp('after obj')    


FO.clickPoint = pointTo;
if FO.EegUpdateBoolean
%     disp('before updateEEG')
    FO = updateEEG(FO,pointTo);
%     disp('after updateEEG')    
end


%Set GUI
% disp('before FO = UpdateText(FO);')
FO = UpdateGUI(FO);
guidata(obj, FO);
% dbstack
% disp('END')
end

function MouseScroll(obj, src)
st = dbstack;%make sure not overloading thestateeditor
if length(st)>1
    for idx = 2:length(st)
        if strmatch (st(idx).file,'TheStateEditor.m')%if state editor is already running
            disp('You scrolled while TheStateEditor was still running, click later')
            return
        end
    end
end

% obj = findobj('tag','StateEditorMaster');  
FO = guidata(obj);

c = get(FO.fig, 'CurrentPoint');
xG = (c(1) >= FO.xplotLims(1)) & (c(1) <=  FO.xplotLims(2));
yG = sum((c(2) >= FO.yplotLims(:, 1)) & c(2) <= FO.yplotLims(:, 2)) == 1;
scrollLoc = 0;
if ~(xG && yG)
    yG = sum((c(2) >= FO.yplotLFPLims(:, 1)) & c(2) <= FO.yplotLFPLims(:, 2)) ~= 0;
    if ~(xG && yG)
        return;
    else
        scrollLoc = 1;
    end
end
switch scrollLoc
    case 0
        d = (c(1) - FO.xplotLims(1))./diff(FO.xplotLims);
        xl = get(FO.sax{1}, 'XLim');
        pointTo = xl(1) + (diff(xl)*d);
        
        if src.VerticalScrollCount > 0
            xl = xl*1.25;
            X1 = ((pointTo) - diff(xl)/2);
            X2 = ((pointTo) + diff(xl)/2);
            if X1 < FO.lims(1)
                X2 = X2 + (FO.lims(1) - X1);
                X1 = FO.lims(1);
            end
            
            if X2 > FO.lims(2)
                X1 = X1 - (X2 - FO.lims(2));
                X2 = FO.lims(2);
            end
            if X1 < FO.lims(1)
                X1 = FO.lims(1);
            end
            set(FO.sax{1}, 'XLim', [X1, X2]);
        else
            xl = xl*0.75;
            X1 = ((pointTo) - diff(xl)/2);
            X2 = ((pointTo) + diff(xl)/2);
            if X1 < FO.lims(1)
                X2 = X2 + (FO.lims(1) - X1);
                X1 = FO.lims(1);
            end
            if X2 > FO.lims(2)
                X1 = X1 - (X2 - FO.lims(2));
                X2 = FO.lims(2);
            end
            set(FO.sax{1}, 'XLim', [X1, X2]);
        end
    case 1
        if src.VerticalScrollCount > 0
            f = find(histc(FO.eegDisplaySeconds, [0, 0.24, 1.99, 4.99, 14.99, 29.9,  60]) == 1);
            switch f
                case 2
                    delta = 0.25;
                case 3
                    delta = 0.5;
                case 4
                    delta = 1;
                case 5
                    delta = 2.5;
                case 6
                    delta = 5;
                case 7
                    delta = 0;
            end
            FO.eegDisplaySeconds = FO.eegDisplaySeconds + delta;
%             guidata(FO.fig, FO);
            
            if FO.EegUpdateBoolean
                FO = updateEEG(FO,mean(get(FO.sax{1}, 'XLim')));
            end
            FO = UpdateGUI(FO);
            return;
        else
            f = find(histc(FO.eegDisplaySeconds, [0, 0.26, 2.1, 5.1, 15.1, 30.1,  60.1]) == 1);
            switch f
                case 1
                    delta = 0;
                case 2
                    delta = 0.25;
                case 3
                    delta = 0.5;
                case 4
                    delta = 1;
                case 5
                    delta = 2.5;
                case 6
                    delta = 5;
            end
            
            FO.eegDisplaySeconds = FO.eegDisplaySeconds - delta;
%             guidata(FO.fig, FO); 
            if FO.EegUpdateBoolean
                FO = updateEEG(FO,mean(get(FO.sax{1}, 'XLim')));
            end
            FO = UpdateGUI(FO);
            return;
        end
end


% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
FO.clickPoint = pointTo;

if FO.EegUpdateBoolean
    FO = updateEEG(FO,pointTo);
end
FO = UpdateGUI(FO);
guidata(FO.fig, FO); 

end



function unMouseClick(e, src)
global isClicking
isClicking = 0;
end

function Nothing(e, src)
a = 0;
end

function FO = addStateLine(FO,location)
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj);
if isempty(FO.startLine) 
    ax = FO.lax;
    yl = get(ax, 'YLim');
    axes(ax);
    hold on;
    FO.startLine{end + 1} = plot([location, location], [yl(1), yl(2)], '-k');
    
    for i = 1:FO.nCh
        ax = FO.sax{i};
        yl = get(ax, 'YLim');
        axes(ax);
        hold on;
        FO.startLine{end + 1} = plot([location, location], [yl(1), yl(2)], '-k');
    end
    
    ax = FO.max;
    yl = get(ax, 'YLim');
    axes(ax);
    hold on;
    FO.startLine{end + 1} = plot([location, location], [yl(1), yl(2)], '-k');
    
    for i = 1:FO.nCh
        ax = FO.eax{i};
        yl = get(ax, 'YLim');
        axes(ax);
        hold on;
        FO.startLine{end + 1} = plot([location, location], [yl(1), yl(2)], '-k');
    end
    
    FO.startLocation = location;
%     guidata(FO.fig, FO); 
    
else    
    
    for i = 1:length(FO.startLine)%delete graphical lines
        delete(FO.startLine{i});
    end    
    FO = addState(location);

    FO.currAction = 'Browse';
    FO.startLocation = [];
    FO.startLine = {};
    set(gcf, 'Name', ['States: ', FO.baseName, ' - Default']);
    
%     guidata(FO.fig, FO); 

    
%     FO.currAction = 'Browse';
%     if ~isempty(FO.startLine)
%         for i = 1:length(FO.startLine)
%             delete(FO.startLine{i});
%         end
%         FO.startLine = {};
% 
%     end
%     set(gcf, 'Name', ['States: ', FO.baseName, ' - Default']);

%     obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
end
% FO = UpdateText(FO);
guidata(FO.fig, FO); 
end

function FO = addState(loc2)
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
s = FO.currentState;
f1 = dsearchn(FO.to, FO.startLocation);
f2 = dsearchn(FO.to, loc2);

% f = [min([f1, f2]), max([f1, f2])];
f = sort([f1 f2]);
newState = zeros(1, diff(f) + 1) + s;
oldState = FO.States(f(1):f(2));
if length(FO.stateHistory) > FO.stateHistoryNum
    FO.stateHistory = FO.stateHistory(1:FO.stateHistoryNum);
    FO.newStates = FO.newStates(1:FO.stateHistoryNum);
    b = msgbox('Losing a bit of history');
    uiwait(b);
end
% if FO.startLocation > loc2
%     FO.Transitions = [FO.Transitions; s, loc2, FO.startLocation];
% else
    FO.Transitions = [FO.Transitions; s, f(1), f(2)];
% end
FO.TransHistoryTracker = [FO.TransHistoryTracker, 1];
FO.stateHistory{end + 1}.location = f(1);
FO.stateHistory{end}.state = oldState;
FO.stateHistoryNum = FO.stateHistoryNum + 1;
FO.States(f(1):f(2)) = newState;
FO.newStates{end + 1}.state = newState;
FO.newStates{end}.location = f(1);
FO.startLocation  = [];
guidata(FO.fig, FO); 
FO = modifyStates(f(1), newState);
FO = guidata(FO.fig);
% updateEEG;

%Reset StateAxis ative click toggle based on checkbox status - since adding
%state temporarily turns it off.
% FO = ResetStateAxisClicabilityByCheckboxStatus(FO);
FO = UpdateGUI(FO);

end


function FO = modifyStates(startLoc, newState, varargin)

if isempty(varargin)
    makeGrey = 1;
else
    makeGrey = varargin{1};
end
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 

colors = FO.colors;
rows = FO.rows;
newC = ones(size(FO.SM, 1), length(newState), 3);
if makeGrey == 1
    newC = repmat(colors.grey, [size(FO.SM, 1), length(newState), 1]);
end

for i = 1:5
    f = find(newState == i);
    newC(rows{i}, f, :) = repmat(colors.states{i}, [length(rows{i}), length(f), 1]);
end

loc = startLoc:(startLoc + length(newState) - 1);

FO.SM(:, loc, :) = newC;
FO.madeChanges = 1;
set(FO.ilab, 'CData', FO.SM);

FO = ResetStateAxisClicabilityByCheckboxStatus(FO);

guidata(FO.fig, FO); 
end


function undoChange(varargin)
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 

if FO.stateHistoryNum >= 1
    newState = FO.stateHistory{FO.stateHistoryNum}.state;
    startLoc = FO.stateHistory{FO.stateHistoryNum}.location;
    FO.stateHistoryNum = FO.stateHistoryNum - 1;
    guidata(FO.fig, FO); 
    modifyStates(startLoc, newState);
    obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
    latestTrans = max(find(FO.TransHistoryTracker == 1));
    FO.TransHistoryTracker(latestTrans) = 0;
    set(gcf, 'Name', ['States: ', FO.baseName, '- History rewound to ', int2str(FO.stateHistoryNum), ' of ', int2str(length(FO.stateHistory))]);
    
else
    b = msgbox('Sorry pal, we are at the beginning of time, there are no changes to undo');
    uiwait(b);
end
guidata(FO.fig, FO); 
end

function redoChange(varargin)
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj);
if length(FO.newStates) > FO.stateHistoryNum
    FO.stateHistoryNum = FO.stateHistoryNum + 1;
    newState = FO.newStates{FO.stateHistoryNum}.state;
    startLoc = FO.newStates{FO.stateHistoryNum}.location;
    
    guidata(FO.fig, FO); 
    modifyStates(startLoc, newState);
    
    latestTrans = max(find(FO.TransHistoryTracker == 0));
    FO.TransHistoryTracker(latestTrans) = 1;
    guidata(FO.fig, FO); 
    set(gcf, 'Name', ['States: ', FO.baseName, '- History moved forward to ', int2str(FO.stateHistoryNum), ' of ', int2str(length(FO.stateHistory))]);
    
else
    b = msgbox('The management regretfully informs you there are no changes to redo.');
    uiwait(b);
end

guidata(FO.fig, FO); 
end


function changeEEGDisplaySeconds(obj,~)
FO = guidata(obj);
sec = get(FO.eegDisplayWidthBox,'String');
sec = str2double(sec);
if sec>30
    disp('Cannot exceed 30 seconds width')
    sec = 30;
end
set(FO.eegDisplayWidthBox,'String',num2str(sec))
FO.eegDisplaySeconds = sec;
FO = updateEEG(FO);

guidata(FO.fig, FO); 
end



function saveRawEEG(obj) %save eeg data to a file so you don't have to have the .lfp file
% obj = findobj('tag','StateEditorMaster');  
FO = guidata(obj);

RawEegData.data = getappdata(obj,'eeg');
for eix = 1:length(RawEegData.data)
    RawEegData.data{eix} = int16(RawEegData.data{eix});
end

RawEegData.channels = FO.Chs;
save(fullfile(FO.basePath,[FO.baseName, '.RawEEG.eegstates.mat']), 'RawEegData');
disp(['EEG/LFP saved to disk as ' FO.baseName, '.RawEEG.eegstates.mat'])
end

function ToggleUpdateEeg(obj,~)
% obj = findobj('tag','StateEditorMaster');  
FO = guidata(obj);

for i = 1:FO.nCh%blank out plots
    set(FO.Eplot{i}, 'XData', []);
end
FO.EegUpdateBoolean = ~FO.EegUpdateBoolean;
guidata(FO.fig, FO); 
end

function ToggleStateAxisActive(obj,~)
% obj = findobj('tag','StateEditorMaster');  
FO = guidata(obj);
FO.stateAxisToggleModeBool = ~FO.stateAxisToggleModeBool;
guidata(FO.fig, FO); 
end


function saved = saveStates
% ?Give user question of whether to impose minimum MA onto
% SleepState.state.mat?

%get figure and figure data
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj);
FO = guidata(obj(end));
baseName = FO.baseName;
basePath = FO.basePath;

%Load baseName.SleepState.states.mat
[SleepState,sleepstatefilename] = bz_LoadStates(basePath,'SleepState');

%Check if the SleepState file is new, manually detected, or auto-detected
if isempty(SleepState)
    STATESFILETYPE = 'new';
    SleepState.detectorinfo.detectorname = 'TheStateEditor';
    SleepState.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd'); 
    SleepState.idx.statenames = {'','','','',''};
    SleepState.idx.timestamps = FO.to;
elseif isfield(FO,'AutoScore')
    STATESFILETYPE = 'auto';
elseif isfield(SleepState,'detectorinfo')
        if isfield(SleepState.detectorinfo,'detectorname')
            switch SleepState.detectorinfo.detectorname
                case 'SleepScoreMaster'
                    STATESFILETYPE = 'auto';
                case 'TheStateEditor'
                    STATESFILETYPE = 'manual';
                otherwise
                    STATESFILETYPE = 'unknown'; 
            end
        end
else
	STATESFILETYPE = 'unknown';   
end


%Put the new states in format for buzcode
%Make an buzcode-style idx structure
if isfield(SleepState,'idx')
    idx.statenames = SleepState.idx.statenames;
else
    idx.statenames = {'WAKE','','NREM','','REM'}; %assume...
end
%Interpolate back to the scoring timestamps
idx.timestamps = SleepState.idx.timestamps;
idx.states = interp1(FO.to,FO.States,idx.timestamps,'nearest');


%Make a buzcode-style ints structure (should wrap this into IDXtoINT.mat)
ints = bz_IDXtoINT(idx,'nameStates',true);


switch STATESFILETYPE
    case 'auto'
        %Save old Autoscoring in AutoScoreInts
        if ~isfield(SleepState,'AutoScoreInts') 
            if isfield(SleepState,'detectorinfo')
                if isfield(SleepState.detectorinfo,'detectorname')
                    if strcmp(SleepState.detectorinfo.detectorname,'SleepScoreMaster')
                        disp('Original State Scoring from SleepScoreMaster detected...')
                        disp('   saving old states as SleepState.AutoScoreInts')
                        SleepState.AutoScoreInts = SleepState.ints;
                    end
                end
            end
        end

        %Save histsandthreshs... different depending on whether AutoScore was
        %used or not
        HistAndThreshAlready_Bool = 0;
        if isfield(FO,'AutoScore')
            if isfield(FO.AutoScore,'histsandthreshs')
                HistAndThreshAlready_Bool = 1;
            end
        end
        if HistAndThreshAlready_Bool
            histsandthreshs = FO.AutoScore.histsandthreshs;
        else
            answer = questdlg({'We usually save hisograms and thresholds for sleep autoscoring, but they are not here... Regenerate them?'});
            switch answer
                case 'Yes'
                    histsandthreshs = SSHistogramsAndThresholds_In(baseName,basePath);
            end
        end
        if exist('histsandthreshs','var')
            SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs = histsandthreshs;
        end
end

%Write the new ints/idx
SleepState.ints = ints;
SleepState.idx = idx;
SleepState.detectorinfo.LastManualUpdate = datestr(now,'yyyy-mm-dd');

%Save the results!
save(sleepstatefilename,'SleepState')


%Make a new figure
try
    ClusterStates_MakeFigure(SleepState,basePath,true);
    disp('Figures Saved to StateScoreFigures')
catch
    disp('Figure making error')
end

%If autoscored, calculate and save new StateEpisodes
switch STATESFILETYPE
    case 'auto'
        StatesToEpisodes(SleepState,basePath);
end



b = msgbox(['Saved work to ', baseName, '.SleepState.states.mat']);
saved = 1;
uiwait(b);

% 
% global answer1
% saved = 0;
% name = [];
% name = FO.baseName;
% name = [name '-states'];
% 
% answer1 = 0;
% FO.saveFig = figure('Position', [382   353   438   200]);
% tx1 = annotation('textbox',  'Position', [0.02, 0.8, 0.5, 0.1], 'string', 'Please Enter File Name:', 'EdgeColor', 'none');
% name1 = uicontrol('style', 'edit', 'string', name, 'FontSize', 10, 'Units', 'Normalized', 'Position', [0.37, 0.78, 0.5, 0.12]);
% set(name1, 'HorizontalAlignment', 'left');
% 
% incEvents = uicontrol('style', 'checkbox', 'string', 'Include Event Times');
% set(incEvents, 'Units', 'normalized', 'Position', [0.1, 0.6, 0.35, 0.1], 'Value', 1);
% 
% incTransitions = uicontrol('style', 'checkbox', 'string', 'Include Transition Times (higher resolution than state vector)');
% set(incTransitions, 'Units', 'normalized', 'Position', [0.1, 0.48, 0.72, 0.1], 'Value', 1);
% 
% incHist = uicontrol('style', 'checkbox', 'string', 'Include History of Changes');
% set(incHist, 'Units', 'normalized', 'Position', [0.1, 0.36, 0.35, 0.1], 'Value', 1);
% 
% saveb = uicontrol('style', 'pushbutton', 'string', 'Save', 'Callback', 'global answer1; uiresume(gcbf); answer1 = 1;');
% set(saveb, 'Units', 'normalized', 'Position', [0.4, 0.1, 0.25, 0.2], 'FontSize', 12);
% 
% cancelb = uicontrol('style', 'pushbutton', 'string', 'Cancel', 'Callback', 'global answer1; uiresume(gcbf); answer1 = 0;');
% set(cancelb, 'Units', 'normalized', 'Position', [0.7, 0.1, 0.25, 0.2], 'FontSize', 12);
% 
% uiwait(FO.saveFig);
% fileName = get(name1, 'string');
% includeH = get(incHist, 'Value');
% includeEvents = get(incEvents, 'Value');
% incTransitions = get(incTransitions, 'Value');
% 
% close(FO.saveFig);
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
% if answer1 == 0
%     return;
% else
%     oldFile = 0;
%     try
%         p = fopen([fileName, '.mat']);
%         fclose(p);
%         oldFile = 1;
%     catch
%         oldFie = 0;
%         
%     end
%     if oldFile == 1
%         choice = questdlg([fileName, '.mat already exists, do you wish to overwrite?']);
%         if strcmpi(choice, 'Yes');
%             delete([fileName, '.mat']);
%         else
%             return;
%         end
%     end
%     
%     states = FO.States;
%     events = FO.Events;
%     
%     toSave = [' ''states'','];
%     if includeEvents == 1
%         toSave = [toSave, ' ''events'','];
%     end
%     
%     if includeEvents == 1
%         if ~isempty(FO.TransHistoryTracker)
%             transitions = FO.Transitions(FO.TransHistoryTracker == 1, :);
%         else
%             transitions = [];
%         end
%         toSave = [toSave, ' ''transitions'','];
%     end
%     
%     if includeH
%         history.stateHistory = FO.stateHistory;
%         history.newStates = FO.newStates;
%         history.stateHistoryNum = FO.stateHistoryNum;
%         toSave = [toSave, ' ''history'','];
%     end
%     
%     try
%         eval(['save(fileName, ', toSave(1:(end - 1)), ');']);
%         b = msgbox(['Saved work to ', fileName, '.mat']);
%         saved = 1;
%         uiwait(b);
%         FO.madeChanges = 0;
%         guidata(FO.fig, FO); 
%     catch
%         b = msgbox(['Warning, failed to save ', fileName, '.mat']);
%         saved = 0;
%         uiwait(b);
%         
%     end
%     
%     
% end
% 
% guidata(FO.fig, FO); 
end


function LoadStates(filepath)
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
global answer1;
answer1 = 0;
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); 
% if answer1 == 0
%     return;
% else
    
    if exist('filepath','var')
        [path,name] = fileparts(filepath);
    else
        [name, path] = uigetfile('*mat', 'Choose a file to load:');
    end
    
    if name == 0
        return;
    end
    
    if strcmp(name(end-21:end),'.SleepState.states.mat')%buzcode format
        thispath = FO.basePath;
        timevector = FO.to;
        states = bz_LoadStates_StateEditorWrapper_In(thispath,timevector);
        FO.States = states;
        
    elseif strcmp(name(end-14:end),'_SleepScore.mat')%2016 format
        load([path,name])
%         stateslen = max([max(max(StateIntervals.NREMstate)) max(max(StateIntervals.REMstate)) max(max(StateIntervals.WAKEstate)) ]); 
        stateslen = size(FO.to,1);
        states = zeros(1,stateslen);
        states(find(inttoboolIn(SleepState.ints.WAKEstate))) = 1;
        %states(find(inttoboolIn(SleepState.ints.MAstate))) = 2;
        states(find(inttoboolIn(SleepState.ints.NREMstate))) = 3;
        states(find(inttoboolIn(SleepState.ints.REMstate))) = 5;
        states = cat(2,states,zeros(1,numel(FO.States)-length(states)));
        FO.States = states;
    else %Andres original TheStateEditor format

        FO.loadFig = figure('Position', [382   353   438   200], 'Name', 'Load');
        warn1 = {'\color{red}\fontsize{12}Warning: \color{black}Loading files will overwrite current work.'};
        tx1 = annotation('textbox',  'Position', [0.02, 0.9, 0.9, 0.1], 'string', warn1 , 'EdgeColor', 'none');

        loadStates = uicontrol('style', 'checkbox', 'string', 'Load State Vector (if ''.states'' field exists)');
        set(loadStates, 'Units', 'normalized', 'Position', [0.1, 0.71, 0.72, 0.15], 'Value', 1, 'FontSize', 10);

        loadEvents = uicontrol('style', 'checkbox', 'string', 'Load Event Matrix (if ''.events'' field exists)');
        set(loadEvents, 'Units', 'normalized', 'Position', [0.1, 0.51, 0.72, 0.15], 'Value', 1, 'FontSize', 10);

        loadTransitions =  uicontrol('style', 'checkbox', 'string', 'Load Transition Matrix (if ''.transitions'' field exists)');
        set(loadTransitions, 'Units', 'normalized', 'Position', [0.1, 0.31, 0.72, 0.15], 'Value', 1, 'FontSize', 10);

        loadb = uicontrol('style', 'pushbutton', 'string', 'Load ''.mat'' File', 'Callback', 'global answer1; uiresume(gcbf); answer1 = 1;');
        set(loadb, 'Units', 'normalized', 'Position', [0.25, 0.06, 0.4, 0.2], 'FontSize', 12);

        cancelb = uicontrol('style', 'pushbutton', 'string', 'Cancel', 'Callback', 'global answer1; uiresume(gcbf); answer1 = 0;');
        set(cancelb, 'Units', 'normalized', 'Position', [0.7, 0.06, 0.25, 0.2], 'FontSize', 12);
        uiwait(FO.loadFig);

        loadStates = get(loadStates, 'Value');
        loadEvents = get(loadEvents, 'Value');
        loadTransitions = get(loadTransitions, 'Value');

        close(FO.loadFig);

        
        newS = load([path, name]);

        if ~isstruct(newS)
            warndlg('Input must be a structure with fields ''.states'', ''.events'' and/or ''.transitions''.')
            return;
        end

        loaded = {};
        if isfield(newS, 'States')
            st = 'States';
        else
            st = 'states';
        end

        if loadStates == 1    
            if isfield(newS, st)
                if length(size(newS.states))==2 && sum(size(newS.states)==1) && numel(newS.states) == numel(FO.States)%make tolerant to vert or horiz vectors
                    newS.states = newS.states(:)';
                end

                if sum(size(FO.States) == size(newS.(st))) ~= 2
                    b = msgbox({'Error: states field must be a 1xN vector file', 'where N == the number of bins.'});
                    uiwait(b);
                    obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
                    return;
                end

                FO.States = newS.(st);
                loaded{end + 1} = 'states vector';
            end
        end

        if loadEvents == 1    
            if isfield(newS, 'events')
                events = newS.events;
                if size(events, 2) ~= 2 & ~isempty(events)
                    b = msgbox({'Error: events field must be a Nx2 matrix file', 'where N == the number of events.'});
                    uiwait(b);
                    obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
                    return;
                end

                FO.Events = events;
                loaded{end + 1} = 'event matrix';
            end
        end


        if loadTransitions == 1
            if isfield(newS, 'transitions')
                transitions = newS.transitions;
                if size(transitions, 2) ~= 3 & ~isempty(transitions)
                    b = msgbox({'Error: transitions field must be a Nx3 matrix file', 'where N == the number of transitions.'});
                    uiwait(b);
                    obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
                    return;
                end

                FO.Transitions = transitions;
                loaded{end + 1} = 'transition matrix';
            end
        end
    end
% end

% successMsg = ['Found and loaded '];
% 
% for i = 1:length(loaded)
%     successMsg = [successMsg, loaded{i}, ', '];
% end
% successMsg = [successMsg(1:(end - 2)), '.'];
% 
% guidata(FO.fig, FO); 
% b = msgbox(successMsg);
% uiwait(b);

EN = FO.eventNum;
if isempty(FO.Events)
    FO = updateEventLines(FO,[]);
else
    FO = updateEventLines(FO,FO.Events(FO.Events(:, 1) == EN, 2));
end
modifyStates(1, FO.States, 0);
guidata(obj, FO); 
end

% 
% function LoadStatesAutoNoMsgs
% % global answer1;
% % answer1 = 0;
% % FO.loadFig = figure('Position', [382   353   438   200], 'Name', 'Load');
% % warn1 = {'\color{red}\fontsize{12}Warning: \color{black}Loading files will overwrite current work.'};
% % tx1 = annotation('textbox',  'Position', [0.02, 0.9, 0.9, 0.1], 'string', warn1 , 'EdgeColor', 'none');
% 
% % loadStates = uicontrol('style', 'checkbox', 'string', 'Load State Vector (if ''.states'' field exists)');
% % set(loadStates, 'Units', 'normalized', 'Position', [0.1, 0.71, 0.72, 0.15], 'Value', 1, 'FontSize', 10);
% % 
% % loadEvents = uicontrol('style', 'checkbox', 'string', 'Load Event Matrix (if ''.events'' field exists)');
% % set(loadEvents, 'Units', 'normalized', 'Position', [0.1, 0.51, 0.72, 0.15], 'Value', 1, 'FontSize', 10);
% % 
% % loadTransitions =  uicontrol('style', 'checkbox', 'string', 'Load Transition Matrix (if ''.transitions'' field exists)');
% % set(loadTransitions, 'Units', 'normalized', 'Position', [0.1, 0.31, 0.72, 0.15], 'Value', 1, 'FontSize', 10);
% % 
% % loadb = uicontrol('style', 'pushbutton', 'string', 'Load ''.mat'' File', 'Callback', 'global answer1; uiresume(gcbf); answer1 = 1;');
% % set(loadb, 'Units', 'normalized', 'Position', [0.25, 0.06, 0.4, 0.2], 'FontSize', 12);
% % 
% % cancelb = uicontrol('style', 'pushbutton', 'string', 'Cancel', 'Callback', 'global answer1; uiresume(gcbf); answer1 = 0;');
% % set(cancelb, 'Units', 'normalized', 'Position', [0.7, 0.06, 0.25, 0.2], 'FontSize', 12);
% % uiwait(FO.loadFig);
% % 
% % loadStates = get(loadStates, 'Value');
% % loadEvents = get(loadEvents, 'Value');
% % loadTransitions = get(loadTransitions, 'Value');
% loadStates = 1;
% loadEvents = 1;
% loadTransitions = 1;
% 
% % close(FO.loadFig);
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
% % if answer1 == 0
% %     return;
% % else
% %     
% %     [name, path] = uigetfile('*mat', 'Choose a state vector to load:');
% %     
% %     if name == 0
% %         return;
% %     end
%     
%     path = cd;
%     name = [FO.baseName '-states.mat'];
%     newS = load(fullfile(path, name));
%     
%     if ~isstruct(newS)
%         warndlg('Input must be a structure with fields ''.states'', ''.events'' and/or ''.transitions''.')
%         return;
%     end
%     
%     loaded = {};
%     if isfield(newS, 'States')
%         st = 'States';
%     else
%         st = 'states';
%     end
%     
%     if loadStates == 1    
%         if isfield(newS, st)
% 
%             if sum(size(FO.States) == size(newS.(st))) ~= 2
%                 b = msgbox({'Error: states field must be a 1xN vector file', 'where N == the number of bins.'});
%                 uiwait(b);
%                 obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
%                 return;
%             end
% 
%             FO.States = newS.(st);
% %             loaded{end + 1} = 'states vector';
%         end
%     end
%     
%     if loadEvents == 1    
%         if isfield(newS, 'events')
%             events = newS.events;
%             if size(events, 2) ~= 2 & ~isempty(events)
%                 b = msgbox({'Error: events field must be a Nx2 matrix file', 'where N == the number of events.'});
%                 uiwait(b);
%                 obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
%                 return;
%             end
% 
%             FO.Events = events;
% %             loaded{end + 1} = 'event matrix';
%         end
%     end
%     
%     
%     if loadTransitions == 1
%         if isfield(newS, 'transitions')
%             transitions = newS.transitions;
%             if size(transitions, 2) ~= 3 & ~isempty(transitions)
%                 b = msgbox({'Error: transitions field must be a Nx3 matrix file', 'where N == the number of transitions.'});
%                 uiwait(b);
%                 obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
%                 return;
%             end
%             
%             FO.Transitions = transitions;
% %             loaded{end + 1} = 'transition matrix';
%         end
%     end
% % end
% 
% % successMsg = ['Found and loaded '];
% % 
% % for i = 1:length(loaded)
% %     successMsg = [successMsg, loaded{i}, ', '];
% % end
% % successMsg = [successMsg(1:(end - 2)), '.'];
% 
% guidata(FO.fig, FO); 
% % b = msgbox(successMsg);
% % uiwait(b);
% 
% EN = FO.eventNum;
% if isempty(FO.Events)
%     FO = updateEventLines(FO,[]);
% else
%     FO = updateEventLines(FO,FO.Events(FO.Events(:, 1) == EN, 2));
% end
% if isfield(newS,'states')
%     modifyStates(1, newS.(st), 0);
% end
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
% FO.madeChanges = 0;
% guidata(FO.fig, FO); 
% end



function out = convWithIn(sp, win1)
%function out = convWith(sp, win1) - win1 convtrimmed with collumns of sp

out = [];
for I = 1:size(sp, 2)
    
    out = [out, convtrimIn(sp(:, I), win1)/sum(win1)];
end
end

function FO = UpdateGUI(FO);
if ~exist('FO','var')
    obj = findobj('tag','StateEditorMaster');  FO = guidata(obj);
end
action = FO.currAction;

%reset action of clicking state axes plot
set(FO.actionDisp, 'color', 'k','backgroundcolor','none');%defaults to offset ResetState

switch action
    case 'Browse'
        set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', ' ', '\fontsize{20}Browse'});
        set(FO.fig,'Pointer','hand');
    case 'Add'
        set(gcf,'Pointer','arrow');
        colors = FO.colors;
        s1 = FO.currentState;
        if s1 == 0
            s2 = 1;
        else
            s2 = s1;
        end
        h = ['\fontsize{42}\color[rgb]{', num2str(colors.states{s2}), '}',int2str(s1)];
        set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', 'Add state:', h});
    case 'AddEvent'
        h = ['\fontsize{42}', int2str(FO.eventNum)];
        set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', 'Add event:', h});
    case 'DeleteEvent'
        h = ['\fontsize{42}', int2str(FO.eventNum)];
        set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', 'Delete event:', h});
    case 'Zoom'
        set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', ' ', '\fontsize{20}Zooming', '\fontsize{20}about'});
        set(gcf,'Pointer','cross');
    case 'FreqResize'
        set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', ' ', '\fontsize{20}Rescale Frequencies'});
        set(gcf,'Pointer','hand');
    case 'ResetState'
        set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', ' ', '\fontsize{20}Choose New State (1-5)'});
        set(FO.actionDisp, 'color', [1 0 0],'backgroundcolor',[0 0 1]);
        set(gcf,'Pointer','arrow');
end

set(FO.max,'xticklabel',num2str(get(FO.max,'xtick')'));

drawnow
set(FO.lastClickDisp, 'String', {'Last Click at sec:', num2str(FO.clickPoint, 7), ['(of ', num2str(FO.lims(2), 7), ')']});
set(FO.eegWidthDisp, 'String', ['\bf\color{red}\fontsize{11}', num2str(FO.eegDisplaySeconds), ' sec']);
% if isempty(FO.startLocation)
%     set(FO.startLocDisp, 'Visible', 'off');
% else
% %     set(FO.startLocDisp, 'String', {'First bound at sec:', num2str(FO.startLocDisp, 7)});
%     set(FO.startLocDisp, 'Visible', 'on');
% end

set(FO.xlimbox, 'String', int2str(round(diff(get(FO.sax{1}, 'XLim')))));
% guidata(FO.fig, FO); 
figure(FO.fig)
end

function ResizeFreqY(direction)
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj);
spec = getappdata(gcf,'spec');

if direction == 1
    m = FO.maxFreq + 10;
    if m <= max(FO.fo);
        FO.maxFreq = m;
        for i = 1:FO.nCh
            set(FO.iSpec{i}, 'CData', spec{i}(:, FO.fo <= m)', 'YData', FO.fo(FO.fo <= m));
            set(FO.sax{i}, 'Ylim', [min(FO.fo), m]);
        end
    end
else
    m = FO.maxFreq - 10;
    if m >= 10
        FO.maxFreq = m;
        for i = 1:FO.nCh
            set(FO.iSpec{i}, 'CData', spec{i}(:, FO.fo <= m)', 'YData', FO.fo(FO.fo <= m));
            set(FO.sax{i}, 'Ylim', [min(FO.fo), m]);
        end
    end
end

guidata(FO.fig, FO); 
end

function CloseDialog(e, src)
FO = guidata(e);
if isempty(FO)
    delete(gcf);
    return;
end
if FO.madeChanges == 0
    warning('on', 'MATLAB:hg:default_child_strategy:IllegalPermutation')
    delete(e);
else
    choice = questdlg(['Save all your hard work before closing?']);
    switch choice
        case 'Yes'
            saved = saveStates;
            if saved == 1
                warning('on', 'MATLAB:hg:default_child_strategy:IllegalPermutation')
                delete(e);
            else
                return;
            end
        case 'No'
            warning('on', 'MATLAB:hg:default_child_strategy:IllegalPermutation')
            delete(e);
        case 'Cancel'
            return;
    end
end

end


function ChangeSmoothingWindow(e, src)
FO = guidata(e);
spec = getappdata(gcf,'spec');
unsmoothedSpec = getappdata(gcf,'unsmoothedSpec');

val = get(FO.hanningWDisp, 'Value');

Woptions = FO.Woptions;
newW = Woptions(val);
if FO.hanningW == newW
    return;
else
    FO.hanningW = newW;
    if newW == 0
        for i = 1:FO.nCh
            spec{i} = log10(unsmoothedSpec{i});
        end
    else
        for i = 1:FO.nCh
            spec{i} = log10(convWithIn(unsmoothedSpec{i}, hanning(FO.hanningW)));
        end
    end
    for i = 1:FO.nCh
        set(FO.iSpec{i}, 'CData', spec{i}(:, FO.fo <= FO.maxFreq)');
    end
end

setappdata(gcf,'spec',spec)
guidata(FO.fig, FO); 
end

function OverlayDisplay(e, src)
FO = guidata(e);
unsmoothedSpec = getappdata(gcf,'unsmoothedSpec');

switch get(FO.overlayDisp, 'Value')
    case 1  %none
        if isempty(FO.overlayLines)
            return;
        else
            for i = 1:length(FO.overlayLines)
                delete(FO.overlayLines{i})
            end
            FO.overlayLines = {};
        end
    case 2  %Theta ratio
        if ~isempty(FO.overlayLines)
            for i = 1:length(FO.overlayLines)
                delete(FO.overlayLines{i})
            end
            FO.overlayLines = {};
        end
        
        maxF = FO.maxFreq;
        fo = FO.fo;
        for i = 1:length(FO.sax)
            m = mean(unsmoothedSpec{i}(:, fo >= 5 & fo <= 10), 2)./mean(unsmoothedSpec{i}(:, fo >= 0.5 & fo <= 4), 2);
            if FO.hanningW > 0
                m = convtrimIn(m, hanning(FO.hanningW));
            end
            m = m - prctile(m, 1);
            m = m./prctile(m, 99);
            range = maxF*(1/2);
            base = maxF*(1/2);
            m  = m*range;
            m = m + base;
            axes(FO.sax{i});
            hold on;
            FO.overlayLines{i} = plot(FO.to, m, '-w', 'LineWidth', 2.5);
        end
    case 3  %From SleepScoreMaster
        basePath = FO.basePath;
        SleepState = bz_LoadStates(basePath,'SleepState');   
        
        maxF = FO.maxFreq;
        
        if isempty(SleepState.detectorinfo.detectionparms.SleepScoreMetrics)
            disp('No SleepScoreMetrics to Overlay')
        else
            broadbandSlowWave = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.broadbandSlowWave;
            thratio = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;
            over_timestamps = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
            chans = FO.Chs;

            overlaychoicefig = figure('closerequestfcn',@OverlaySleepStateSelectCallback);
            for cidx = 1:length(chans)
                   bg(cidx) = uibuttongroup(overlaychoicefig,...
                      'Position',[(cidx-1)*.33 0 .3 1],...
                      'Title',['Overlay for Ch' num2str(chans(cidx))]);

                    % Create radio buttons in the button group.
                    r1(cidx) = uicontrol(bg(cidx),'Style','radiobutton',...
                          'String','Broaband SlowWave Power',...
                          'Units','Normalized',....
                          'Position',[.05 .55 1 .1]);
                    r2(cidx) = uicontrol(bg(cidx),'Style','radiobutton',...
                          'String','Theta Ratio 5-10Hz/2-20Hz',...
                          'Units','Normalized',....
                          'Position',[.05 .1 1 .08]);

                    if chans(cidx) == SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID
                        r1(cidx).Value = true;
                        r2(cidx).Value = false;
                    end
                    if chans(cidx) == SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID
                        r1(cidx).Value = false;
                        r2(cidx).Value = true;
                    end

            end            
            closebutt = uicontrol('style','pushbutton','units','normalized',... %lol butt. - Good Dan  - :D 
                'position',[.91 .05 .08 .1],'String','Finish',...
                'callback',@OverlaySleepStateSelectCallback);
            localguidata = v2struct(overlaychoicefig,bg,r1,r2,closebutt);
            guidata(overlaychoicefig,localguidata)
            waitfor(overlaychoicefig)

            choices = get(FO.fig,'userdata');
            for cidx = 1:length(choices)%for each channel/choice (should be same)
                switch choices(cidx)
                    case 1
                        t = broadbandSlowWave;
                    case 2
                        t = thratio;
                end

                t = t';
                t = t - prctile(t, 1);
                t = t./prctile(t, 99);
                range = maxF*(1/2);
                base = maxF*(1/2);
                t  = t*range;
                t = t + base;

                axes(FO.sax{cidx});
                hold on;
                FO.overlayLines{cidx} = plot(over_timestamps, t, '-w', 'LineWidth', 2.5);
            end
        end                        
  
        
    case 4  %From File
        helpdlg({['load a .mat with a single variable with n columns of time bins'],...
            ['(n = ', int2str(length(FO.to)),') and up to ', int2str(FO.nCh), ' rows. Successive rows of the'],...
            ['input will be displayed overlayed on on successive'],...
            ['spectrogram channels']});
        
        [name, path] = uigetfile('*mat', 'Choose overlay data to load:');
        
        maxF = FO.maxFreq;
        

        if name == 0
            guidata(FO.fig, FO); 
            set(FO.overlayDisp, 'Value', 1);
        else

            input1 = load([path, name]);

            if isstruct(input1) && ~isfield(input1,'timestamps')
                t = fieldnames(input1);
                input1 = input1.(t{1});
            end
            
            %Buzcode structure with timestamps and data - still needs
            %work....
            if isfield(input1,'timestamps')
                maxF = FO.maxFreq;
                
                [ m ] = bz_NormToRange(input1.data,[0.5*maxF maxF]);
                hold on;
                FO.overlayLines{1} = plot(input1.timestamps, m, '-w', 'LineWidth', 2.5);
            else


                if size(input1, 2) ~= length(FO.to)
                    b = msgbox('Error: number of columns in input does not match the number of bins');
                    uiwait(b);

                    set(FO.overlayDisp, 'Value', 1);
                    guidata(FO.fig, FO); 
                    return;
                end

                if ~isempty(FO.overlayLines)
                    for i = 1:length(FO.overlayLines)
                        delete(FO.overlayLines{i})
                    end
                    FO.overlayLines = {};
                end

                m1 = min([FO.nCh; size(input1, 1)]);
                maxF = FO.maxFreq;
                for i = 1:m1
                    m = input1(i, :);
                    m = m - prctile(m, 1);
                    m = m./prctile(m, 99);
                    range = maxF*(1/2);
                    base = maxF*(1/2);
                    m  = m*range;
                    m = m + base;
                    axes(FO.sax{i});
                    hold on;
                    FO.overlayLines{i} = plot(FO.to, m, '-w', 'LineWidth', 2.5);
                end
            
            end
        end
                
end

guidata(FO.fig, FO); 
FO = UpdateGUI(FO);

end

function OverlaySleepStateSelectCallback(obj,ev)
f = get(obj,'parent');
lgd = guidata(f);

for idx = 1:length(lgd.r1) %only BBSlowWave and ThetaRatio
    if lgd.r1(idx).Value
        out(idx) = 1;%if value1 is 1, save out as 1
    else
        out(idx) = 2;%else save as 2
    end
end
TSEFig = findobj('tag','StateEditorMaster');  
set(TSEFig,'userdata',out)
delete(f)

end


function motion = LoadFromWhl(baseName, tos, varargin)
if length(varargin) >= 1
    whl = varargin{1};
else
    [whl,to,~] = LoadFromWhlHelper1(baseName);
end
if size(whl,2) == 2 %if two LEDs then take the one with fewest NaN's, then fill in NaN's from other light if neccesary
    d = [abs(diff(whl(:, 1))) + abs(diff(whl(:, 2))), abs(diff(whl(:, 3))) + abs(diff(whl(:, 4)))];
    [~, mGood] = find(min(mean(isnan(d))));
    if mGood == 1
        mBad = 2;
        
    else
        mBad = 1;
    end
    motion1 = d(:, mGood);
    motion1(isnan(d(:, mGood))) = d(isnan(d(:, mGood)), mBad);
else
    motion1 = abs(diff(whl));
end
f = find(isnan(motion1));
motion = motion1;

for i = 1:length(f) %%%this looks 3 seconds before and after each NaN value to find a non-NaN estimate (mean of non-NaN neighbors)
    h = [f(i) - 39*3, f(i) + 39*3];
    if h(1) < 1
        h(1) = 1;
    else
        if h(2) > length(motion1)
            h(2) = length(motion1);
        end
    end
    c = motion1(h(1):h(2));
    motion(f(i)) = mean(c(~isnan(c)));
end

to2 = to(2:end) - (diff(to)/2);

motion2 = motion;
motion = [];
for i = 1:length(tos) %avearge over same bins as spectrogram
    motion = [motion, mean(motion2(to2 > tos(i) & to2 <= (tos(i) + 1)))];
end


end

function [whl,t,GoodRanges] = LoadFromWhlHelper1(fbasename)
% USAGE
% [whl,t,GoodRanges,ep] = LoadPosition(fbasename)
%
% output:
%   whl: the position matrix
%   t: time vector
%   Good Ranges:

Fs = 1250/32;

whlt = dlmread([fbasename '.whl']);
[whl GoodRanges] = LoadFromWhlHelper2(whlt);

t = (1:size(whlt,1))'/Fs;
GoodRanges = GoodRanges/Fs;
end

function  [cWhl, GoodRanges_F] = LoadFromWhlHelper2(Whl, StretchLen, JumpSize, Gap)

% If the gap between the good strech is more than StrethcLen in terms of Whl row number,remove interporated values.
if nargin<2
    StretchLen = 30;
end

% If the Gap between the good strech is more than JumpSize, remove interporated values.
if nargin<3
    JumpSize = 30;
end

% if the distance between the two contimous rows are more than Gap centimeter, It's a big jump and do not use as an input for inpterp1.
if nargin<4,
    Gap = 30;
end

nWhl = size(Whl,1);

% interpolate missing values or large jumps.
% the value of whl(:,3:4) is also taken into account for the rang of interpolation.
% A transision to and form (-1,-1) should be taken as a BigJump.

% I hsould use distance, not the one dimentinal projection of trajectory, by the way.

whltemp = Whl;
whltemp(find(whltemp)==-1) = -Gap;
dist_F = sqrt(diff(whltemp(:,1)).^2+diff(whltemp(:,2)).^2);
dist_R = sqrt(diff(whltemp(:,3)).^2+diff(whltemp(:,4)).^2);
BigJump_F = dist_F>Gap;
BigJump_R = dist_R>Gap;

Good_F = find(Whl(:,1)>-1 & ~([BigJump_F;0] | [0;BigJump_F]));
Bad_F = find(~(Whl(:,1)>-1 & ~([BigJump_F;0] | [0;BigJump_F])));
Good_R = find(Whl(:,3)>-1 & ~([BigJump_R;0] | [0;BigJump_R]));
Bad_R = find(~(Whl(:,3)>-1 & ~([BigJump_R;0] | [0;BigJump_R])));

whltemp(Bad_F,1:2) = -Gap;
whltemp(Bad_R,3:4) = -Gap;

WhlNaN = Whl;
WhlNaN(find(Whl==-1)) = NaN;

% Give -1 outside of the interpolation.

if length(Good_F)<2 || length(Good_R)<2;
    cWhl(:,1:2) = -ones(size(Whl,1),2);
else
    cWhl(:,1:2) = interp1(Good_F, Whl(Good_F,1:2), 1:nWhl, 'linear', -1);
    cWhl(:,3:4) = interp1(Good_R, Whl(Good_R,3:4), 1:nWhl, 'linear', -1);
end

% find missing stretches for Front LED
dGoodF = [-(whltemp(1,1)==-Gap) ; diff(whltemp(:,1)>-Gap)];
BadStartF = find(dGoodF<0);
BadEndF = find(dGoodF>0)-1;
% if last point is bad, need to finish off BadEnd
if Whl(end,1)==-1
    BadEndF = [BadEndF; nWhl];
end

if length(BadStartF)>length(BadEndF)
    BadEndF = [BadEndF; nWhl];
end


% find ranges to chuck
% jump size ...
if any(BadStartF>0)
    
    StartIndF = clip(BadStartF-1, 1, nWhl); % StartInd and EndInd give the
    EndIndF = clip(BadEndF+1, 1, nWhl);     % points you are interpolating between
    
    dist_F = sqrt((Whl(StartIndF,1)-Whl(EndIndF,1)).^2+(Whl(StartIndF,2)-Whl(EndIndF,2)).^2);
    ToChuckF = find(BadEndF-BadStartF>=StretchLen ...
        | dist_F > JumpSize);
    % chuck em
    
    for i=ToChuckF(:)'
        cWhl(BadStartF(i):BadEndF(i),1:2) = NaN;
    end
end

% find missing stretches for Rear LED
dGoodR = [-(whltemp(1,3)==-Gap) ; diff(whltemp(:,3)>-Gap)];
BadStartR = find(dGoodR<0);
BadEndR = find(dGoodR>0)-1;
% if last point is bad, need to finish off BadEnd
if Whl(end,3)==-1
    BadEndR = [BadEndR; nWhl];
end

if length(BadStartR)>length(BadEndR)
    BadEndR = [BadEndR; nWhl];
end


% find ranges to chuck
% jump size ...
if any(BadStartR>0)
    StartIndR = clip(BadStartR-1, 1, nWhl); % StartInd and EndInd give the
    EndIndR = clip(BadEndR+1, 1, nWhl);     % points you are interpolating between
    
    dist_R = sqrt((Whl(StartIndR,3)-Whl(EndIndR,3)).^2+(Whl(StartIndR,4)-Whl(EndIndR,4)).^2);
    ToChuckR = find(BadEndR-BadStartR>=StretchLen ...
        | dist_R > JumpSize);
    
    % chuck em
    for i=ToChuckR(:)'
        cWhl(BadStartR(i):BadEndR(i),3:4) = NaN;
    end
end


if 0 % OLD VERSION (BUG?)
    % % now find good ranges
    % dcGood = [-(Whl(1,1)==1) ; diff(cWhl(:,1)>-1)];
    % GoodStart = find(dcGood>0);
    % GoodEnd = find(dcGood<0)-1;
    % % if last point is good, need to finish GoodEnd
    % if cWhl(end,1)>-1
    %     GoodEnd = [GoodEnd; nWhl];
    % end
    % GoodRanges = [GoodStart, GoodEnd];
else
    dcGood_F = diff([0; cWhl(:,1)>-1; 0]);
    GoodStart_F = find(dcGood_F>0);
    GoodEnd_F = find(dcGood_F<0)-1;
    GoodRanges_F = [GoodStart_F, GoodEnd_F];
end


return

end


function motion = LoadTimeStampValuePairs(tos,fname,varname);
if ~exist('varname','var')
    t = load(fname);
else
    t = load(fname,varname);
end
fn = fieldnames(t);
t = getfield(t,fn{1});
vals = t(:,2);
times = t(:,1);
motion = ResampleTolerant_IN(vals,length(tos),length(times));

end


function Par = LoadParIn(FileName, varargin)
% LoadPar(FileName)
% loads the specified par file and returns a structure with these elements:
%
% .FileName      -> name of file loaded from
% .nChannels     -> number of total channels
% .nBits         -> number of bits of the file
% .SampleTime    -> time, in microseconds, of 1 sample (ie 1e6 / sample rate)
% .HiPassFreq    -> High pass filter frequency
% .nElecGps      -> number of electrodes (i.e. electrode groups)
% .ElecGp        -> a cell array giving the channels in the electrodes
%                    e.g. if .ElectrodeGroup{3} = [2 3 4 5], electrode 3
%                    is a tetrode for channels 2 3 4 and 5. 
% channel numbers here are from 0. be carefull.
[SpecInfo] = DefaultArgsIn(varargin,{1});

if ~isempty(strfind(FileName,'.par'))
    FileBase = FileName(1:strfind(FileName,'.par')-1);
elseif ~isempty(strfind(FileName,'.xml'))
    FileBase = FileName(1:strfind(FileName,'.xml')-1);
else 
    FileBase = FileName;
end


if FileExistsIn([FileBase '.xml']) %& ~isempty(strfind(FileName,'.xml'))
    Par = LoadParameters([FileBase '.xml']);
elseif FileExistsIn([FileBase '.par'])


    % open file

    fp = fopen([FileBase '.par'], 'r');
    Par.FileName = FileBase;

    % read in nChannels and nBits
    Line = fgets(fp);
    A = sscanf(Line, '%d %d');
    Par.nChannels = A(1);
    Par.nBits = A(2);

    % read in SampleTime and HiPassFreq
    Line = fgets(fp);
    A = sscanf(Line, '%d %f', 2);
    Par.SampleTime = A(1);
    Par.HiPassFreq = A(2);

    % read in nElectrodes
    Line = fgets(fp);
    if Line==-1
        fclose(fp);
        return;
    end
    A = sscanf(Line, '%d', 1);
    Par.nElecGps = A(1);

    % read in ElectrodeGroup
    for i=1:Par.nElecGps
        Line = fgets(fp);
        A = sscanf(Line, '%d');
        Par.ElecGp{i} = A(2:end);
    end
    fclose(fp);
else
    error('Par or Xml file do not exist!');
end

if SpecInfo
    if FileExistsIn([FileBase '.eeg.par'])
        if ~isfield(Par,'nElecGps')
            ParTmp = LoadPar([FileBase '.par']);
        else
            ParTmp = Par;
        end

        EegPar=LoadEegPar(FileBase);
        for el=1:ParTmp.nElecGps
            for eegel=1:EegPar.nElec
                if ~isempty(intersect(ParTmp.ElecGp{el},EegPar.ElecChannels{eegel}))
                    Par.ElecLoc{el} = EegPar.ElecLoc{eegel};
                end
            end
        end
   
    end
end
end


function [c] = convtrimIn(a,b)
% CONVTRIM trimmed convolution
% c = convtrim(a,b) convolves vectors A and B. The resulting
% vector is length LENGTH(a)
%
% this function is a wrapper for conv - the only difference is the trimming

if (length(a) <= length(b))
    error('convtrim: the length of vector a must be larger than vector b');
end

tempC = conv(a,b);
FrontTrim = floor(length(b)/2);

if (mod(length(b),2) ~= 0)
    BackTrim = floor(length(b)/2);
else
    BackTrim = floor(length(b)/2)-1;
end

c = tempC(FrontTrim+1:end-BackTrim);



end

function changeXlim(obj,~)
FO = guidata(obj);
xmax = diff(FO.lims);
n1 = get(FO.xlimbox, 'String');
n2 = round(str2double(n1));
oldx = get(FO.sax{1}, 'XLim');
if (n2 <= xmax) && n2 > 0
    m1 = mean(oldx);
    newx = [m1 - n2/2, m1 + n2/2];
    if min(newx) < 0
        newx = newx - min(newx);
    end
    if max(newx) > FO.lims(2)
        newx = [FO.lims(2) - diff(newx), FO.lims(2)];
    end
    
    set(FO.sax{1}, 'XLim', newx);
    axes(FO.sax{1});
    FO = UpdateGUI(FO);
    if FO.EegUpdateBoolean
        FO = updateEEG(FO);
    end
else
    set(FO.xlimbox, 'String', int2str(round(diff(oldx))));
    return;
end
end


function goToSecond(obj,~)
% get xlims to generate range
% get the string to get what the mean should be
% make newx based on that
% if less than 0, or greater than max possible, keep the range and set the
%    one end at the boundary and the other
% set FO.sax
FO = guidata(obj);%get guidata
% xmax = diff(FO.lims);%max x possible
gotopoint = str2double(get(FO.gotosecondbox, 'String'));
oldrange = get(FO.sax{1}, 'XLim');
windowwidth = diff(oldrange);
% if (n2 <= xmax) & n2 > 0

newrange = [gotopoint - windowwidth/2, gotopoint + windowwidth/2];

if min(newrange) < 0
    newrange = [0 windowwidth-1];
end
if max(newrange) > FO.lims(2)
    newrange = [(FO.lims(2)-windowwidth+1), FO.lims(2)];
end

% set(FO.sax{1}, 'XLim', newrange);  %action step: set axis1, if there are other axes, updateEEG will set any other axes to match it
% axes(FO.sax{1});
for idx = 1:length(FO.sax)
    set(FO.sax{idx}, 'XLim', newrange);  %action step: set axis1, if there are other axes, updateEEG will set any other axes to match it
%     axes(FO.sax{idx});
end

set(FO.gotosecondbox, 'String', '');

FO = UpdateGUI(FO);
if FO.EegUpdateBoolean
    if gotopoint > FO.lims(1) && gotopoint <= FO.lims(2)
        FO = updateEEG(FO,gotopoint);
    end
end
% else
%     set(FO.xlimbox, 'String', int2str(round(diff(oldx))));
%     return;
% end
end


function rightScreenOver(obj,~)
% move over with only 10% overlap from previous screen
FO = guidata(obj);%get guidata
oldrange = get(FO.sax{1}, 'XLim');
windowwidth = diff(oldrange);

newrange = oldrange+windowwidth*0.9;

if min(newrange) < 0
    newrange = [0 windowwidth-1];
end
if max(newrange) > FO.lims(2)
    newrange = [(FO.lims(2)-windowwidth+1), FO.lims(2)];
end

for idx = 1:length(FO.sax)
    set(FO.sax{idx}, 'XLim', newrange);  %action step: set axis1, if there are other axes, updateEEG will set any other axes to match it
end

set(FO.gotosecondbox, 'String', '');

FO = UpdateGUI(FO);
if FO.EegUpdateBoolean
    FO = updateEEG(FO);
end
end

function leftScreenOver(obj,event)
% move over with only 10% overlap from previous screen
FO = guidata(obj);%get guidata
oldrange = get(FO.sax{1}, 'XLim');
windowwidth = diff(oldrange);

newrange = oldrange-windowwidth*0.9;

if min(newrange) < 0
    newrange = [0 windowwidth-1];
end
if max(newrange) > FO.lims(2)
    newrange = [(FO.lims(2)-windowwidth+1), FO.lims(2)];
end

for idx = 1:length(FO.sax)
    set(FO.sax{idx}, 'XLim', newrange);  %action step: set axis1, if there are other axes, updateEEG will set any other axes to match it
end

set(FO.gotosecondbox, 'String', '');

FO = UpdateGUI(FO);
if FO.EegUpdateBoolean
    FO = updateEEG(FO);
end
end

function FO = ChangeThisStateAssigmnent_Setup(FO,pointTo)
states = FO.States;
ds = logical(abs(diff(states)));%find state switches with 1's
ds = cat(2,ds(:)',0);%pad

%find start
spanstart = find(ds(1:pointTo)>0,1,'last')+1;
if isempty(spanstart)
    spanstart = 1;
end
if spanstart<1
    spanstart = 1;
end
%find end
spanend = find(ds(pointTo:end)>0,1,'first') + pointTo-1;
if isempty(spanend)
    spanend = length(states);
end
if spanend>length(states)
    spanend = length(states);
end
span = [spanstart spanend];
FO.spanToReset = span;

FO.currAction = 'ResetState';
% set(FO.actionDisp, 'String', {'\fontsize{12}\bfCurrent Action:', ' ', '\fontsize{20}Choose New State (1-5)'});
set(FO.fig,'KeyReleaseFcn', {@ChangeThisStateAssigmnent_Key});
end

function ChangeThisStateAssigmnent_Key(obj,ev)
FO = guidata(obj);
statetype = ev.Key;
if sum(strcmp(statetype,{'1','2','3','4','5'}))%if key was 1:5
    if ~isempty(FO.spanToReset)
        statetype = str2double(statetype);        
        span = FO.spanToReset;
        newState = statetype * ones(1,diff(span)+1);
        
% s = FO.currentState;
% f1 = dsearchn(FO.to, FO.startLocation);
% f2 = dsearchn(FO.to, pointTo);
% f = [min([f1, f2]), max([f1, f2])];
% newState = zeros(1, diff(f) + 1) + s;

        oldState = FO.States(span(1):span(2));
        if length(FO.stateHistory) > FO.stateHistoryNum
            FO.stateHistory = FO.stateHistory(1:FO.stateHistoryNum);
            FO.newStates = FO.newStates(1:FO.stateHistoryNum);
            b = msgbox('Losing a bit of history');
            uiwait(b);
        end
        FO.States(span(1):span(2)) = newState;

        FO.Transitions = [FO.Transitions; statetype, span(1), span(2)];
        FO.TransHistoryTracker = [FO.TransHistoryTracker, 1];
        FO.stateHistory{end + 1}.location = span(1);
        FO.stateHistory{end}.state = oldState;
        FO.stateHistoryNum = FO.stateHistoryNum + 1;
        FO.newStates{end + 1}.state = newState;
        FO.newStates{end}.location = span(1);
        FO.startLocation  = [];
% % guidata(FO.fig, FO); 
% FO = modifyStates(f(1), newState);
% FO = guidata(FO.fig);
% % updateEEG;
% 
% %Reset StateAxis ative click toggle based on checkbox status - since adding
% %state temporarily turns it off.
% % FO = ResetStateAxisClicabilityByCheckboxStatus(FO);
% FO = UpdateGUI(FO);        
        
        guidata(FO.fig, FO); 
        modifyStates(span(1), newState);
        FO = guidata(FO.fig);
    end    
end

FO.spanToReset = [];
FO.currAction = 'Browse';
set(FO.fig,'KeyReleaseFcn', {@DefKey});

FO = UpdateGUI(FO);
guidata(FO.fig,FO)
end

function FO = ResetStateAxisClicabilityByCheckboxStatus(FO)
%Reset StateAxis ative click toggle based on checkbox status - since adding
%state temporarily turns it off.

% if strmatch(FO.currAction,'Add')
    SAACBool = get(FO.stateAxisActiveCheckbox,'value');
    if SAACBool
        FO.stateAxisToggleModeBool = logical(1);
    elseif ~SAACBool 
        FO.stateAxisToggleModeBool = logical(0);
    end
% end

end

function a = clip(b, min, max)
% clip(b, min, max)
% takes a matrix and replaces any elements below min
% with min and any above max with max

shape = size(b);
b = b(:);
b(find(b<min)) = min;
b(find(b>max)) = max;
a = reshape(b, shape);
end

function EventNumber(e, src)

obj = findobj('tag','StateEditorMaster');  
FO = guidata(obj); ;
EN = get(FO.eventDisp, 'Value') - 1;

oldNum = FO.eventNum;
if EN > 0
    FO.eventNum = EN;
else
    FO.eventNum = 'none';
end
guidata(FO.fig, FO); 

% if ~isempty(FO.Events)
    FO = updateEventLines(FO,FO.Events(FO.Events(:, 1) == EN, 2));
% else 
% end

guidata(FO.fig, FO); 
if strcmp(FO.currAction, 'AddEvent')
    if EN == 0
        FO.currAction = 'Browse';
        set(gcf, 'Pointer', 'hand');
        guidata(FO.fig, FO); 
    end
    FO = UpdateGUI(FO);
end

end

function FO = addEvent(FO,location)
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
EN = FO.eventNum;
FO.Events = [FO.Events; EN, location];
% guidata(FO.fig, FO); 
FO = updateEventLines(FO,FO.Events(FO.Events(:, 1) == EN, 2));
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
% FO.CurrEventLines{end + 1} = [];
% 
% for i = 1:FO.nCh
%     ax = FO.sax{i};
%     yl = [0 300];
%     axes(ax);
%     hold on;
%     FO.CurrEventLines{end}{panelN} = plot([location, location], [yl(1), yl(2)], ':m', 'LineWidth', 2);
%     panelN = panelN + 1;
%     ax = FO.eax{i};
%     yl = [-10000 10000];
%     axes(ax);
%     hold on;
%     FO.CurrEventLines{end}{panelN} = plot([location, location], [yl(1), yl(2)], ':m', 'LineWidth', 2);
%     panelN = panelN + 1;
% end
% 
% ax = FO.max;
% yl = get(FO.max, 'YLim');
% axes(ax);
% hold on;
% FO.CurrEventLines{end}{panelN} = plot([location, location], [yl(1), yl(2)], ':m', 'LineWidth', 2);
% FO.madeChanges = 1;

FO.currAction = 'Browse';
set(FO.fig, 'Pointer', 'hand');
% guidata(FO.fig, FO); 
if FO.EegUpdateBoolean
    FO = updateEEG(FO,location);
end
end

function FO = deleteEvent(FO,location)
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
x1 = get(gca, 'XLim');
EN = FO.eventNum;
events = FO.Events;
ind1 = find(events(:, 1) == EN);
id = dsearchn(events(events(:, 1) == EN, 2), location);
if (abs(location - events(ind1(id), 2))/diff(x1)) < 0.01;
   
    events = events((1:size(events, 1)) ~= ind1(id), :);
    FO.Events = events;
%     guidata(FO.fig, FO); 
    FO = updateEventLines(FO,events(events(:, 1) == EN, 2));
%     obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
    
    FO.currAction = 'Browse';
    set(FO.fig, 'Pointer', 'hand');
end

% guidata(FO.fig, FO); 

if FO.EegUpdateBoolean
    FO = updateEEG(FO,location);
end
% FO = UpdateText(FO);;
end

function FO = updateEventLines(FO,newTimes)
% obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
if isempty(newTimes)
    for i = 1:length(FO.CurrEventLines)
        set(FO.CurrEventLines{i}, 'XData', [], 'YData', []);
    end
    FO.CurrEventLines = {};
else
    panelN = 1;
    
    
    
    z = zeros(length(newTimes)*2, 1);
    z(mod(1:(length(newTimes)*2), 2) == 0) = newTimes;
    z(mod(1:(length(newTimes)*2), 2) == 1) = newTimes;
    if mod(length(newTimes), 2) == 1
        y = repmat([0; 1; 1; 0], (length(newTimes) - 1)/2, 1);
        y = [y; 0; 1];
    else
        y = repmat([0; 1; 1; 0], (length(newTimes))/2, 1);
    end
    
    linesExist = length(FO.CurrEventLines) > 0;
    for i = 1:FO.nCh
        
        axes(FO.sax{i});
        hold on;
        yl = y;
        yl(yl == 0) = FO.saxYLim(1) - 10;
        yl(yl == 1) = FO.saxYLim(2) + 10;
        if linesExist == 0
            FO.CurrEventLines{panelN} = plot(z, yl , ':m', 'LineWidth', 2);
        else
            set(FO.CurrEventLines{panelN}, 'XData', z, 'YData', yl);
        end
        panelN = panelN + 1;
        axes(FO.eax{i});
        hold on;
        yl = y;
        yl(yl == 0) = FO.eegYLim(1) - 10;
        yl(yl == 1) = FO.eegYLim(2) + 10;
        if linesExist == 0
            FO.CurrEventLines{panelN} = plot(z, yl , ':m', 'LineWidth', 2);
        else
            set(FO.CurrEventLines{panelN}, 'XData', z, 'YData', yl);
        end
        panelN = panelN + 1;
    end
    axes(FO.max);
    hold on;
    yl = y;
    yl(yl == 0) = FO.mpYLim(1) - 10;
    yl(yl == 1) = FO.mpYLim(2) + 10;
    if linesExist == 0
        FO.CurrEventLines{panelN} = plot(z, yl , ':m', 'LineWidth', 2);
    else
        set(FO.CurrEventLines{panelN}, 'XData', z, 'YData', yl);
    end
    panelN = panelN + 1;
end
newString = 'none';

events = FO.Events;
if isempty(events)
for i = 1:10
        newString = [newString, '|', int2str(i), ' (0 events)'];
    end
else
    for i = 1:10
        newString = [newString, '|', int2str(i), ' (', int2str(sum(events(:, 1) == i)), ' events)'];
    end
end
set(FO.eventDisp, 'String', newString);


% guidata(FO.fig, FO); 

end

function nextEvent
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
events = FO.Events;
EN = FO.eventNum;
if isempty(events)
    warndlg('No events currently selected.');
    return;
end

ev = events(events(:, 1) == EN, 2);
if length(ev) == 0
    warndlg('No events currently selected.');
    return;
end


ev = sort(ev);
currPoint = mean(get(FO.eax{1}, 'XLim'));
nextE = ev(find(ev > currPoint, 1, 'first'));
if isempty(nextE)
    warndlg('No further events to jump to.');
    return;
end

saxLim = diff(get(FO.sax{1}, 'XLim'));

newLims = [nextE - saxLim/2, nextE + saxLim/2];
if newLims(1) < FO.lims(1)
    newLims(2) = newLims(2) + (FO.lims(1) - newLims(1)); 
    newLims(1) = FO.lims(1);
end

if newLims(2) > FO.lims(2)
    newLims(1) = newLims(1) + (FO.lims(2) - newLims(2));
    newLims(2) = FO.lims(2);
end

if newLims(1) < FO.lims(1)
    newLims(2) = newLims(2) + (FO.lims(1) - newLims(1)); 
    newLims(1) = FO.lims(1);
end

set(FO.sax{1}, 'XLim', newLims);
if FO.EegUpdateBoolean
    FO = updateEEG(FO,nextE);
end
end



function previousEvent
obj = findobj('tag','StateEditorMaster');  FO = guidata(obj); ;
events = FO.Events;
EN = FO.eventNum;
if isempty(events)
    warndlg('No events currently selected.');
    return;
end

ev = events(events(:, 1) == EN, 2);
if length(ev) == 0
    warndlg('No events currently selected.');
    return;
end


ev = sort(ev);
currPoint = mean(get(FO.eax{1}, 'XLim'));
nextE = ev(find(ev < currPoint, 1, 'last'));
if isempty(nextE)
    warndlg('No previous events to jump to.');
    return;
end

saxLim = diff(get(FO.sax{1}, 'XLim'));

newLims = [nextE - saxLim/2, nextE + saxLim/2];
if newLims(1) < FO.lims(1)
    newLims(2) = newLims(2) + (FO.lims(1) - newLims(1)); 
    newLims(1) = FO.lims(1);
end

if newLims(2) > FO.lims(2)
    newLims(1) = newLims(1) + (FO.lims(2) - newLims(2));
    newLims(2) = FO.lims(2);
end

if newLims(1) < FO.lims(1)
    newLims(2) = newLims(2) + (FO.lims(1) - newLims(1)); 
    newLims(1) = FO.lims(1);
end

set(FO.sax{1}, 'XLim', newLims);
if FO.EegUpdateBoolean
    FO = updateEEG(FO,nextE);
end
end

function E = FileExistsIn(name)
if length(dir(name)) > 0
    E = 1;
else
    E = 0;
end
end

function skullCursor

skull = [NaN   NaN   NaN     2     1     1     1     1     1     1     1     2   NaN   NaN   NaN   NaN
    NaN     2     1     1     1     1   NaN     1     1     1     1     1     1     2   NaN   NaN
    2     1     1     1     1   NaN   NaN   NaN     1     1     1     1     1     1     2   NaN
    2     1     1     1     1     1   NaN     1     1     1     1     1     1     1     2   NaN
    2     1     1     1     1     1     1     1     1     1     1     1     1     1     2   NaN
    2     1     1     2     2     2     1     1     1     2     2     2     1     1     2   NaN
    2     1     1     2     2     2     1     1     1     2     2     2     1     1     2   NaN
    2     1     1     1     2     1     1     1     1     1     2     1     1     1     2   NaN
    NaN     2     1     1     1     1     1     2     1     1     1     1     1     2   NaN   NaN
    NaN     2     1     1     1     1     2     2     2     1     1     1     1     2   NaN   NaN
    NaN   NaN     2     1     1     1     1     1     1     1     1     1     2   NaN   NaN   NaN
    NaN   NaN     2     1     2     1     1     1     1     1     2     1     2   NaN   NaN   NaN
    NaN   NaN     2     1     1     2     1     1     1     2     1     1     2   NaN   NaN   NaN
    NaN   NaN   NaN     2     1     1     2     2     2     1     1     2   NaN   NaN   NaN   NaN
    NaN   NaN   NaN   NaN     2     1     1     1     1     1     2   NaN   NaN   NaN   NaN   NaN
    NaN   NaN   NaN   NaN   NaN     2     1     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN];

set(gcf, 'PointerShapeCData', skull );
set(gcf, 'Pointer', 'custom');
end

function [y, f, t, phi, FStats]=mtchglongIn(varargin);
%function [yo, fo, to, phi, FStats]=mtchglong(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
% Multitaper Time-Frequency Cross-Spectrum (cross spectrogram)
% for long files - splits data into blockes to save memory
% function A=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,nTapers)
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%
% output yo is yo(f, t)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, t, Ch1, Ch2)
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere

% Original code by Partha Mitra - modified by Ken Harris 
% and adopted for long files and phase by Anton Sirota
% Also containing elements from specgram.m

% default arguments and that
[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t] = mtparamIn(varargin);

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFFTChunks,nFreqBins, nChannels, nChannels)); % output array
if nargout>3
    phi=complex(zeros(nFFTChunks,nFreqBins, nChannels, nChannels));
end
nFFTChunksall= nFFTChunks;
freemem = FreeMemoryIn;
BlockSize = 2^8;
nBlocks = ceil(nFFTChunksall/BlockSize);
%h = waitbar(0,'Wait..');
for Block=1:nBlocks
    %   waitbar(Block/nBlocks,h);
    minChunk = 1+(Block-1)*BlockSize;
    maxChunk = min(Block*BlockSize,nFFTChunksall);
    nFFTChunks = maxChunk - minChunk+1;
    iChunks = [minChunk:maxChunk];
    Periodogram = complex(zeros(nFreqBins, nTapers, nChannels, nFFTChunks)); % intermediate FFTs
    Temp1 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    Temp2 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    Temp3 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    eJ = complex(zeros(nFreqBins, nFFTChunks));
    tmpy =complex(zeros(nFreqBins,nFFTChunks, nChannels, nChannels));
    % calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
    [Tapers V]=dpss(WinLength,NW,nTapers, 'calc');
    % New super duper vectorized alogirthm
    % compute tapered periodogram with FFT 
    % This involves lots of wrangling with multidimensional arrays.
    
    TaperingArray = repmat(Tapers, [1 1 nChannels]);
    for j=1:nFFTChunks
        jcur = iChunks(j);
        Segment = x((jcur-1)*winstep+[1:WinLength], :);
        if (~isempty(Detrend))
            Segment = detrend(Segment, Detrend);
        end;
        SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
        TaperedSegments = TaperingArray .* SegmentsArray;
        
        fftOut = fft(TaperedSegments,nFFT);
        normfac = sqrt(2/nFFT); %to get back rms of original units
        Periodogram(:,:,:,j) = fftOut(select,:,:)*normfac; %fft(TaperedSegments,nFFT);
        % Periodogram: size  = nFreqBins, nTapers, nChannels, nFFTChunks
    end	
    if nargout>4
        U0 = repmat(sum(Tapers(:,1:2:end)),[nFreqBins,1,nChannels,   nFFTChunks]);
        Mu = sq(sum(Periodogram(:,1:2:end,:,:) .* conj(U0), 2) ./  sum(abs(U0).^2, 2));
        Num = abs(Mu).^2;
        Sp = sq(sum(abs(Periodogram).^2,2));
        chunkFS = (nTapers-1) * Num ./ (Sp ./ sq(sum(abs(U0).^2, 2))- Num );
        %	sum(abs(Periodogram - U0.*repmat(Mu,[1,nTapers,1,1])), 2);
        FStats(iChunks, :, :)  = permute(reshape(chunkFS, [nFreqBins, nChannels, nFFTChunks]),[ 3 1, 2]);
    end
    % Now make cross-products of them to fill cross-spectrum matrix
    for Ch1 = 1:nChannels
        for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
            Temp1 = reshape(Periodogram(:,:,Ch1,:), [nFreqBins,nTapers,nFFTChunks]);
            Temp2 = reshape(Periodogram(:,:,Ch2,:), [nFreqBins,nTapers,nFFTChunks]);
            Temp2 = conj(Temp2);
            Temp3 = Temp1 .* Temp2;
            eJ=sum(Temp3, 2);
            tmpy(:,:, Ch1, Ch2)= eJ/nTapers;
            
            % for off-diagonal elements copy into bottom half of matrix
            if (Ch1 ~= Ch2)
                tmpy(:,:, Ch2, Ch1) = conj(eJ) / nTapers;
            end            
            
        end
    end
    
    for Ch1 = 1:nChannels
        for Ch2 = 1:nChannels % don't compute cross-spectra twice
            
            if (Ch1 == Ch2)
                % for diagonal elements (i.e. power spectra) leave unchanged
                y(iChunks,:,Ch1, Ch2) = permute(tmpy(:,:,Ch1, Ch2),[2 1 3 4]);
            else
                % for off-diagonal elements, scale
                
                y(iChunks,:,Ch1, Ch2) = permute((abs(tmpy(:,:,Ch1, Ch2).^2) ...
                    ./ (tmpy(:,:,Ch1,Ch1) .* tmpy(:,:,Ch2,Ch2))), [2 1 3 4]);
                if nargout>3
                    phi(iChunks,:,Ch1,Ch2) = permute(angle(tmpy(:,:,Ch1, Ch2) ...
                        ./ sqrt(tmpy(:,:,Ch1,Ch1) .* tmpy(:,:,Ch2,Ch2))), [2 1 3 4]); 
                end
            end
        end
    end
    
    
end
%close(h);
% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
    % take abs, and use image to display results
    newplot;
    for Ch1=1:nChannels, for Ch2 = 1:nChannels
            subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
	    if Ch1==Ch2
		if length(t)==1
			imagesc([0 1/f(2)],f,20*log10(abs(y(:,:,Ch1,Ch2))+eps)');axis xy; colormap(jet);
		else
			imagesc(t,f,20*log10(abs(y(:,:,Ch1,Ch2))+eps)');axis xy; colormap(jet);
		end
	    else
	    	imagesc(t,f,(abs(y(:,:,Ch1,Ch2)))');axis xy; colormap(jet);
	    end
        end; end;
    xlabel('Time')
    ylabel('Frequency')
end
end


%LoadBinary - Load data from a binary file.
%
%  USAGE
%
%    data = LoadBinary(filename,<options>)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'duration'    duration to read (in s) (default = Inf)
%     'frequency'   sampling rate (in Hz) (default = 20kHz)
%     'start'       position to start reading (in s) (default = 0)
%     'nChannels'   number of data channels in the file (default = 1)
%     'channels'    channels to read (default = all)
%     'precision'   sample precision (default = 'int16')
%     'skip'        number of bytes to skip after each value is read
%                   (default = 0)
%    =========================================================================

% Copyright (C) 2004-2006 by Micha?l Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function data = LoadBinaryIn(filename,varargin)

% Default values
start = 0;
nChannels = 1;
precision = 'int16';
skip = 0;
duration = Inf;
frequency = 20000;
channels = 1;

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help LoadBinary'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help LoadBinary'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'duration',
      duration = varargin{i+1};
      if ~isa(duration,'numeric') | length(duration) ~= 1 | duration < 0,
        error('Incorrect value for property ''duration'' (type ''help LoadBinary'' for details).');
      end
    case 'frequency',
      frequency = varargin{i+1};
      if ~isa(frequency,'numeric') | length(frequency) ~= 1 | frequency <= 0,
        error('Incorrect value for property ''frequency'' (type ''help LoadBinary'' for details).');
      end
    case 'start',
      start = varargin{i+1};
      if ~isa(start,'numeric') | length(start) ~= 1,
        error('Incorrect value for property ''start'' (type ''help LoadBinary'' for details).');
      end
		if start < 0, start = 0; end
    case 'nchannels',
      nChannels = varargin{i+1};
      if ~((round(channels) == channels & channels > 0)) | length(nChannels) ~= 1,
        error('Incorrect value for property ''nChannels'' (type ''help LoadBinary'' for details).');
      end
    case 'channels',
      channels = varargin{i+1};
      if ~(round(channels) == channels & channels > 0)
        error('Incorrect value for property ''channels'' (type ''help LoadBinary'' for details).');
      end
    case 'precision',
      precision = varargin{i+1};
      if ~isa(precision,'char'),
        error('Incorrect value for property ''precision'' (type ''help LoadBinary'' for details).');
      end
    case 'skip',
      skip = varargin{i+1};
      if ~IsPositiveInteger(skip) | length(skip) ~= 1,
        error('Incorrect value for property ''skip'' (type ''help LoadBinary'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help LoadBinary'' for details).']);
  end
end

sizeInBytes = 0;
switch precision,
	case {'uchar','unsigned char','schar','signed char','int8','integer*1','uint8','integer*1'},
		sizeInBytes = 1;
	case {'int16','integer*2','uint16','integer*2'},
		sizeInBytes = 2;
	case {'int32','integer*4','uint32','integer*4','single','real*4','float32','real*4'},
		sizeInBytes = 4;
	case {'int64','integer*8','uint64','integer*8','double','real*8','float64','real*8'},
		sizeInBytes = 8;
end

f = fopen(filename,'r');

% Position file index for reading
start = floor(start*frequency)*nChannels*sizeInBytes;
status = fseek(f,start,'bof');
if status ~= 0,
    fclose(f);
    error('Could not start reading (possible reasons include trying to read past the end of the file).');
end

% Determine number of samples when duration is 'inf'
if isinf(duration),
	fileStart = ftell(f);
	status = fseek(f,0,'eof');
	if status ~= 0,
		fclose(f);
		error('Error reading the data file (possible reasons include trying to read past the end of the file).');
	end
	fileStop = ftell(f);
	nSamplesPerChannel = (fileStop-fileStart)/nChannels/sizeInBytes;
	duration = nSamplesPerChannel/frequency;
	frewind(f);
	status = fseek(f,start,'bof');
	if status ~= 0,
		fclose(f);
		error('Could not start reading (possible reasons include trying to read past the end of the file).');
	end
else
    nSamplesPerChannel = floor(frequency*duration);
    if nSamplesPerChannel ~= frequency*duration,
        %disp(['Warning: rounding duration (' num2str(duration) ' -> ' num2str(nSamplesPerChannel/frequency) ')']);
        duration = nSamplesPerChannel/frequency;
    end
end

% For large amounts of data, read chunk by chunkChangeThisStateAssigmnent_Key

maxSamplesPerChunk = 100000;
nSamples = nChannels*nSamplesPerChannel;
if nSamples > maxSamplesPerChunk,
	% Determine chunk duration and number of chunks
	nSamplesPerChunk = floor(maxSamplesPerChunk/nChannels)*nChannels;
	durationPerChunk = nSamplesPerChunk/frequency/nChannels;
	nChunks = floor(duration/durationPerChunk);
	% Preallocate memory
	data = zeros(nSamplesPerChannel,length(channels),precision);
	% Read all chunks
	i = 1;
	for j = 1:nChunks,
		d = LoadBinaryChunkIn(f,'frequency',frequency,'nChannels',nChannels,'channels',channels,'duration',durationPerChunk,'skip',skip);
		[m,n] = size(d);
		if m == 0, break; end
		data(i:i+m-1,:) = d;
		i = i+m;
%  		h=waitbar(j/nChunks);
	end
%  	close(h)
	% If the data size is not a multiple of the chunk size, read the remainder
	remainder = duration - nChunks*durationPerChunk;
	if remainder ~= 0,
		d = LoadBinaryChunkIn(f,'frequency',frequency,'nChannels',nChannels,'channels',channels,'duration',remainder,'skip',skip);
		[m,n] = size(d);
		if m ~= 0,
			data(i:i+m-1,:) = d;
		end
	end
else
    if skip ~= 0,
        data = fread(f,[nChannels frequency*duration],precision,skip);
    else
        data = fread(f,[nChannels frequency*duration],precision);
    end
    data=data';
    
    if ~isempty(channels),
        data = data(:,channels);
    end
end
fclose(f);
end

% helper function to do argument defaults etc for mt functions
function [x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t,FreqRange] ...
    = mtparamIn(P)

nargs = length(P);

x = P{1};
if (nargs<2 | isempty(P{2})) nFFT = 1024; else nFFT = P{2}; end;
if (nargs<3 | isempty(P{3})) Fs = 1250; else Fs = P{3}; end;
if (nargs<4 | isempty(P{4})) WinLength = nFFT; else WinLength = P{4}; end;
if (nargs<5 | isempty(P{5})) nOverlap = WinLength/2; else nOverlap = P{5}; end;
if (nargs<6 | isempty(P{6})) NW = 3; else NW = P{6}; end;
if (nargs<7 | isempty(P{7})) Detrend = ''; else Detrend = P{7}; end;
if (nargs<8 | isempty(P{8})) nTapers = 2*NW -1; else nTapers = P{8}; end;
if (nargs<9 | isempty(P{9})) FreqRange = [0 Fs/2]; else FreqRange = P{9}; end
% Now do some compuatations that are common to all spectrogram functions
if size(x,1)<size(x,2)
    x = x';
end
nChannels = size(x, 2);
nSamples = size(x,1);

if length(nOverlap)==1
    winstep = WinLength - nOverlap;
    % calculate number of FFTChunks per channel
    %remChunk = rem(nSamples-Window)
    nFFTChunks = max(1,round(((nSamples-WinLength)/winstep))); %+1  - is it ? but then get some error in the chunking in mtcsd... let's figure it later
    t = winstep*(0:(nFFTChunks-1))'/Fs;
else
    winstep = 0;
    nOverlap = nOverlap(nOverlap>WinLength/2 & nOverlap<nSamples-WinLength/2);
    nFFTChunks = length(nOverlap);
    t = nOverlap(:)/Fs; 
end 
%here is how welch.m of matlab does it:
% LminusOverlap = L-noverlap;
% xStart = 1:LminusOverlap:k*LminusOverlap;
% xEnd   = xStart+L-1;
% welch is doing k = fix((M-noverlap)./(L-noverlap)); why?
% turn this into time, using the sample frequency


% set up f and t arrays
if isreal(x)%~any(any(imag(x)))    % x purely real
	if rem(nFFT,2),    % nfft odd
		select = [1:(nFFT+1)/2];
	else
		select = [1:nFFT/2+1];
	end
	nFreqBins = length(select);
else
	select = 1:nFFT;
end
f = (select - 1)'*Fs/nFFT;
nFreqRanges = size(FreqRange,1);
%if (FreqRange(end)<Fs/2)
if nFreqRanges==1
    select = find(f>FreqRange(1) & f<FreqRange(end));
    f = f(select);
    nFreqBins = length(select);
else
    select=[];
    for i=1:nFreqRanges
        select=cat(1,select,find(f>FreqRange(i,1) & f<FreqRange(i,2)));
    end
    f = f(select);
    nFreqBins = length(select);
end
%end
end

%LoadBinaryChunk - Load data chunck from an open binary file.
%
%  USAGE
%
%    data = LoadBinaryChunk(fid,<options>)
%
%    fid            file id (obtained via fopen)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'duration'    duration to read (in s) (default = 1s)
%     'frequency'   sampling rate (in Hz) (default = 20kHz)
%     'start'       position to start reading (in s) (default = from current
%                   index, allowing to read a file chunck by chunck)
%     'nChannels'   number of data channels in the file (default = 1)
%     'precision'   sample precision (default = 'int16')
%    =========================================================================

% Copyright (C) 2004-2006 by Micha?l Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function data = LoadBinaryChunkIn(fid,varargin)

% Default values
start = 0;
fromCurrentIndex = true;
nChannels = 1;
precision = 'int16';
duration = 1;
frequency = 20000;
channels = [];

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help LoadBinaryChunk'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help LoadBinaryChunk'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'duration',
      duration = varargin{i+1};
      if ~isa(duration,'numeric') | length(duration) ~= 1 | duration < 0,
        error('Incorrect value for property ''duration'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'frequency',
      frequency = varargin{i+1};
      if ~isa(frequency,'numeric') | length(frequency) ~= 1 | frequency <= 0,
        error('Incorrect value for property ''frequency'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'start',
      start = varargin{i+1};
      fromCurrentIndex = false;
      if ~isa(start,'numeric') | length(start) ~= 1,
        error('Incorrect value for property ''start'' (type ''help LoadBinaryChunk'' for details).');
      end
		if start < 0, start = 0; end
    case 'nchannels',
      nChannels = varargin{i+1};
      if ~isa(nChannels,'numeric') | length(nChannels) ~= 1,
        error('Incorrect value for property ''nChannels'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'channels',
      channels = varargin{i+1};
      if ~isa(channels,'numeric'),
        error('Incorrect value for property ''channels'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'precision',
      precision = varargin{i+1};
      if ~isa(precision,'char'),
        error('Incorrect value for property ''precision'' (type ''help LoadBinaryChunk'' for details).');
      end
    case 'skip',
      skip = varargin{i+1};
      if ~isa(skip,'numeric') | length(skip) ~= 1,
        error('Incorrect value for property ''skip'' (type ''help LoadBinaryChunk'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help LoadBinaryChunk'' for details).']);
  end
end

sizeInBytes = 0;
switch precision,
	case {'uchar','unsigned char','schar','signed char','int8','integer*1','uint8','integer*1'},
		sizeInBytes = 1;
	case {'int16','integer*2','uint16','integer*2'},
		sizeInBytes = 2;
	case {'int32','integer*4','uint32','integer*4','single','real*4','float32','real*4'},
		sizeInBytes = 4;
	case {'int64','integer*8','uint64','integer*8','double','real*8','float64','real*8'},
		sizeInBytes = 8;
end

% Position file index for reading
if ~fromCurrentIndex,
	start = floor(start*frequency)*nChannels*sizeInBytes;
	status = fseek(fid,start,'bof');
	if status ~= 0,
		error('Could not start reading (possible reasons include trying to read a closed file or past the end of the file).');
	end
end

% Read data chunck
if skip ~= 0,
	data = fread(fid,[nChannels frequency*duration],precision,skip);
else
	data = fread(fid,[nChannels frequency*duration],precision);
end;
data=data';

% Keep only required channels
if ~isempty(channels) & ~isempty(data),
	data = data(:,channels);
end

end

function [y, A] = WhitenSignalIn(x,varargin)

%artype =2; %Signal processing toolbox
artype =1; %arfit toolbox, (crushes sometimes with old version and single data type)

[window,CommonAR, ARmodel,ArOrder] = DefaultArgsIn(varargin,{[],1,[],1});
ArOrder = ArOrder+1;
Trans = 0;
if size(x,1)<size(x,2)
    x = x';
    Transf =1;
end
[nT nCh]  = size(x);
y = zeros(nT,nCh);
if isempty(window)
    seg = [1 nT];
    nwin=1;
else
    nwin = floor(nT/window)+1;
    seg = repmat([1 window],nwin,1)+repmat([0:nwin-1]'*window,1,2);
    if nwin*window>nT
        seg(end,2) =nT;
    end   
end

for j=1:nwin
    if ~isempty(ARmodel) 
        A = ARmodel;
        for i=1:nCh
            y(seg(j,1):seg(j,2),i) = Filter0In(A, x(seg(j,1):seg(j,2),i));
        end
    else
        if CommonAR % meaning common model for all channels and segments!!! 
            for i=1:nCh
                if  j==1 & i==1
                    switch artype
                        case 1
                            [w Atmp] = arfitIn(x(seg(j,1):seg(j,2),i),ArOrder,ArOrder);
                            A = [1 -Atmp];
                        case 2
                            A = arburg(x(seg(j,1):seg(j,2),i),ArOrder);
                    end
                    ARmodel = A;
                end
                y(seg(j,1):seg(j,2),i) = Filter0In(A, x(seg(j,1):seg(j,2),i));
            end
        else
            for i=1:nCh
                switch artype
                    case 1
                        [w Atmp] = arfitIn(x(seg(j,1):seg(j,2),i),ArOrder,ArOrder);
                        A =[1 -Atmp];
                    case 2
                        A = arburg(x(seg(j,1):seg(j,2),i),ArOrder);
                end
                y(seg(j,1):seg(j,2),i) = Filter0In(A, x(seg(j,1):seg(j,2),i));
            end
        end
    end
end

if Trans
    y =y';
end

end


function varargout = DefaultArgsIn(Args, DefArgs)
% auxillary function to replace argument check in the beginning and def. args assigment
% sets the absent or empty values of the Args (cell array, usually varargin)
% to their default values from the cell array DefArgs. 
% Output should contain the actuall names of arguments that you use in the function

% e.g. : in function MyFunction(somearguments , varargin)
% calling [SampleRate, BinSize] = DefaultArgs(varargin, {20000, 20});
% will assign the defualt values to SampleRate and BinSize arguments if they
% are empty or absent in the varargin cell list 
% (not passed to a function or passed empty)
if isempty(Args)
    Args ={[]};
end

if iscell(Args{1}) & length(Args)==1
    Args = Args{1};
end

nDefArgs = length(DefArgs);
nInArgs = length(Args);
%out = cell(nDefArgs,1);
if (nargout~=nDefArgs)
    error('number of defaults is different from assigned');
    keyboard
end
for i=1:nDefArgs
    if (i>nInArgs | isempty(Args{i}))
        varargout(i) = {DefArgs{i}};
    else 
        varargout(i) = {Args{i}};
    end
end

end


function [w, A, C, sbc, fpe, th]=arfitIn(v, pmin, pmax, selector, no_const)
%ARFIT	Stepwise least squares estimation of multivariate AR model.
%
%  [w,A,C,SBC,FPE,th]=ARFIT(v,pmin,pmax) produces estimates of the
%  parameters of a multivariate AR model of order p,
%
%      v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + noise(C),
%
%  where p lies between pmin and pmax and is chosen as the optimizer
%  of Schwarz's Bayesian Criterion. The input matrix v must contain
%  the time series data, with columns of v representing variables
%  and rows of v representing observations.  ARFIT returns least
%  squares estimates of the intercept vector w, of the coefficient
%  matrices A1,...,Ap (as A=[A1 ... Ap]), and of the noise covariance
%  matrix C.
%
%  As order selection criteria, ARFIT computes approximations to
%  Schwarz's Bayesian Criterion and to the logarithm of Akaike's Final
%  Prediction Error. The order selection criteria for models of order
%  pmin:pmax are returned as the vectors SBC and FPE.
%
%  The matrix th contains information needed for the computation of
%  confidence intervals. ARMODE and ARCONF require th as input
%  arguments.
%       
%  If the optional argument SELECTOR is included in the function call,
%  as in ARFIT(v,pmin,pmax,SELECTOR), SELECTOR is used as the order
%  selection criterion in determining the optimum model order. The
%  three letter string SELECTOR must have one of the two values 'sbc'
%  or 'fpe'. (By default, Schwarz's criterion SBC is used.) If the
%  bounds pmin and pmax coincide, the order of the estimated model
%  is p=pmin=pmax. 
%
%  If the function call contains the optional argument 'zero' as the
%  fourth or fifth argument, a model of the form
%
%         v(k,:)' = A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + noise(C) 
%
%  is fitted to the time series data. That is, the intercept vector w
%  is taken to be zero, which amounts to assuming that the AR(p)
%  process has zero mean.

%  Modified 14-Oct-00
%  Authors: Tapio Schneider
%           tapio@gps.caltech.edu
%
%           Arnold Neumaier
%           neum@cma.univie.ac.at

  % n: number of observations; m: dimension of state vectors
  [n,m]   = size(v);     

  if (pmin ~= round(pmin) | pmax ~= round(pmax))
    error('Order must be integer.');
  end
  if (pmax < pmin)
    error('PMAX must be greater than or equal to PMIN.')
  end

  % set defaults and check for optional arguments
  if (nargin == 3)              % no optional arguments => set default values
    mcor       = 1;               % fit intercept vector
    selector   = 'sbc';	          % use SBC as order selection criterion
  elseif (nargin == 4)          % one optional argument
    if strcmp(selector, 'zero')
      mcor     = 0;               % no intercept vector to be fitted
      selector = 'sbc';	          % default order selection 
    else
      mcor     = 1; 		  % fit intercept vector
    end
  elseif (nargin == 5)          % two optional arguments
    if strcmp(no_const, 'zero')
      mcor     = 0;               % no intercept vector to be fitted
    else
      error(['Bad argument. Usage: ', ...
	     '[w,A,C,SBC,FPE,th]=AR(v,pmin,pmax,SELECTOR,''zero'')'])
    end
  end

  ne  	= n-pmax;               % number of block equations of size m
  npmax	= m*pmax+mcor;          % maximum number of parameter vectors of length m

  if (ne <= npmax)
    error('Time series too short.')
  end

  % compute QR factorization for model of order pmax
  [R, scale]   = arqrIn(v, pmax, mcor);

  % compute approximate order selection criteria for models 
  % of order pmin:pmax
  [sbc, fpe]   = arordIn(R, m, mcor, ne, pmin, pmax);

  % get index iopt of order that minimizes the order selection 
  % criterion specified by the variable selector
  [val, iopt]  = min(eval(selector)); 

  % select order of model
  popt         = pmin + iopt-1; % estimated optimum order 
  np           = m*popt + mcor; % number of parameter vectors of length m

  % decompose R for the optimal model order popt according to 
  %
  %   | R11  R12 |
  % R=|          |
  %   | 0    R22 |
  %
  R11   = R(1:np, 1:np);
  R12   = R(1:np, npmax+1:npmax+m);    
  R22   = R(np+1:npmax+m, npmax+1:npmax+m);

  % get augmented parameter matrix Aaug=[w A] if mcor=1 and Aaug=A if mcor=0
  if (np > 0)   
    if (mcor == 1)
      % improve condition of R11 by re-scaling first column
      con 	= max(scale(2:npmax+m)) / scale(1); 
      R11(:,1)	= R11(:,1)*con; 
    end;
    Aaug = (R11\R12)';
    
    %  return coefficient matrix A and intercept vector w separately
    if (mcor == 1)
      % intercept vector w is first column of Aaug, rest of Aaug is 
      % coefficient matrix A
      w = Aaug(:,1)*con;        % undo condition-improving scaling
      A = Aaug(:,2:np);
    else
      % return an intercept vector of zeros 
      w = zeros(m,1);
      A = Aaug;
    end
  else
    % no parameters have been estimated 
    % => return only covariance matrix estimate and order selection 
    % criteria for ``zeroth order model''  
    w   = zeros(m,1);
    A   = [];
  end
  
  % return covariance matrix
  dof   = ne-np;                % number of block degrees of freedom
  C     = R22'*R22./dof;        % bias-corrected estimate of covariance matrix
  
  % for later computation of confidence intervals return in th: 
  % (i)  the inverse of U=R11'*R11, which appears in the asymptotic 
  %      covariance matrix of the least squares estimator
  % (ii) the number of degrees of freedom of the residual covariance matrix 
  invR11 = inv(R11);
  if (mcor == 1)
    % undo condition improving scaling
    invR11(1, :) = invR11(1, :) * con;
  end
  Uinv   = invR11*invR11';
  th     = [dof zeros(1,size(Uinv,2)-1); Uinv];
end

function [R, scale]=arqrIn(v, p, mcor)
%ARQR	QR factorization for least squares estimation of AR model.
%
%  [R, SCALE]=ARQR(v,p,mcor) computes the QR factorization needed in
%  the least squares estimation of parameters of an AR(p) model. If
%  the input flag mcor equals one, a vector of intercept terms is
%  being fitted. If mcor equals zero, the process v is assumed to have
%  mean zero. The output argument R is the upper triangular matrix
%  appearing in the QR factorization of the AR model, and SCALE is a
%  vector of scaling factors used to regularize the QR factorization.
%
%  ARQR is called by ARFIT. 
%
%  See also ARFIT.

%  Modified 29-Dec-99
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu

  % n: number of time steps; m: dimension of state vectors
  [n,m] = size(v);     

  ne    = n-p;                  % number of block equations of size m
  np    = m*p+mcor;             % number of parameter vectors of size m

  % If the intercept vector w is to be fitted, least squares (LS)
  % estimation proceeds by solving the normal equations for the linear
  % regression model
  %
  %                  v(k,:)' = Aaug*u(k,:)' + noise(C)
  %
  % with Aaug=[w A] and `predictors' 
  %
  %              u(k,:) = [1 v(k-1,:) ...  v(k-p,:)]. 
  %
  % If the process mean is taken to be zero, the augmented coefficient
  % matrix is Aaug=A, and the regression model
  %
  %                u(k,:) = [v(k-1,:) ...  v(k-p,:)]
  %
  % is fitted. 
  % The number np is the dimension of the `predictors' u(k). 

  % Assemble the data matrix K (of which a QR factorization will be computed)
  K = zeros(ne,np+m);                 % initialize K
  if (mcor == 1)
    % first column of K consists of ones for estimation of intercept vector w
    K(:,1) = ones(ne,1);
  end
  
  % Assemble `predictors' u in K 
  for j=1:p
    K(:, mcor+m*(j-1)+1:mcor+m*j) = [v(p-j+1:n-j, :)];
  end
  % Add `observations' v (left hand side of regression model) to K
  K(:,np+1:np+m) = [v(p+1:n, :)];
  
  % Compute regularized QR factorization of K: The regularization
  % parameter delta is chosen according to Higham's (1996) Theorem
  % 10.7 on the stability of a Cholesky factorization. Replace the
  % regularization parameter delta below by a parameter that depends
  % on the observational error if the observational error dominates
  % the rounding error (cf. Neumaier, A. and T. Schneider, 2001:
  % "Estimation of parameters and eigenmodes of multivariate
  % autoregressive models", ACM Trans. Math. Softw., 27, 27--57.).
  q     = np + m;             % number of columns of K
  delta = (q^2 + q + 1)*eps;  % Higham's choice for a Cholesky factorization
  scale = sqrt(delta)*sqrt(sum(K.^2));   
  R     = triu(qr([K; diag(scale)]));

% Add `observations' v (left hand side of regression model) to K
K(:,np+1:np+m) = [v(p+1:n, :)];

% Compute regularized QR factorization of K: The regularization
% parameter delta is chosen according to Higham's (1996) Theorem
% 10.7 on the stability of a Cholesky factorization. Replace the
% regularization parameter delta below by a parameter that depends
% on the observational error if the observational error dominates
% the rounding error (cf. Neumaier, A. and T. Schneider, 2001:
% "Estimation of parameters and eigenmodes of multivariate
% autoregressive models", ACM Trans. Math. Softw., 27, 27--57.).
q     = np + m;             % number of columns of K
delta = (q^2 + q + 1)*eps;  % Higham's choice for a Cholesky factorization
scale = sqrt(delta)*sqrt(sum(K.^2));
R     = triu(qr([K; diag(scale)]));
end

function [sbc, fpe, logdp, np] = arordIn(R, m, mcor, ne, pmin, pmax)
%ARORD	Evaluates criteria for selecting the order of an AR model.
%
%  [SBC,FPE]=ARORD(R,m,mcor,ne,pmin,pmax) returns approximate values
%  of the order selection criteria SBC and FPE for models of order
%  pmin:pmax. The input matrix R is the upper triangular factor in the
%  QR factorization of the AR model; m is the dimension of the state
%  vectors; the flag mcor indicates whether or not an intercept vector
%  is being fitted; and ne is the number of block equations of size m
%  used in the estimation. The returned values of the order selection
%  criteria are approximate in that in evaluating a selection
%  criterion for an AR model of order p < pmax, pmax-p initial values
%  of the given time series are ignored.
%
%  ARORD is called by ARFIT. 
%	
%  See also ARFIT, ARQR.

%  For testing purposes, ARORD also returns the vectors logdp and np,
%  containing the logarithms of the determinants of the (scaled)
%  covariance matrix estimates and the number of parameter vectors at
%  each order pmin:pmax.

%  Modified 17-Dec-99
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu

imax 	  = pmax-pmin+1;        % maximum index of output vectors

% initialize output vectors
sbc     = zeros(1, imax);     % Schwarz's Bayesian Criterion
fpe     = zeros(1, imax);     % log of Akaike's Final Prediction Error
logdp   = zeros(1, imax);     % determinant of (scaled) covariance matrix
np      = zeros(1, imax);     % number of parameter vectors of length m
np(imax)= m*pmax+mcor;

% Get lower right triangle R22 of R:
%
%   | R11  R12 |
% R=|          |
%   | 0    R22 |
%
R22     = R(np(imax)+1 : np(imax)+m, np(imax)+1 : np(imax)+m);

% From R22, get inverse of residual cross-product matrix for model
% of order pmax
invR22  = inv(R22);
Mp      = invR22*invR22';

% For order selection, get determinant of residual cross-product matrix
%       logdp = log det(residual cross-product matrix)
logdp(imax) = 2.*log(abs(prod(diag(R22))));

% Compute approximate order selection criteria for models of
% order pmin:pmax
i = imax;
for p = pmax:-1:pmin
    np(i)      = m*p + mcor;	% number of parameter vectors of length m
    if p < pmax
        % Downdate determinant of residual cross-product matrix
        % Rp: Part of R to be added to Cholesky factor of covariance matrix
        Rp       = R(np(i)+1:np(i)+m, np(imax)+1:np(imax)+m);
        
        % Get Mp, the downdated inverse of the residual cross-product
        % matrix, using the Woodbury formula
        L        = chol(eye(m) + Rp*Mp*Rp')';
        N        = L \ Rp*Mp;
        Mp       = Mp - N'*N;
        
        % Get downdated logarithm of determinant
        logdp(i) = logdp(i+1) + 2.* log(abs(prod(diag(L))));
    end
    
    % Schwarz's Bayesian Criterion
    sbc(i) = logdp(i)/m - log(ne) * (ne-np(i))/ne;
    
    % logarithm of Akaike's Final Prediction Error
    fpe(i) = logdp(i)/m - log(ne*(ne-np(i))/(ne+np(i)));
    
    % Modified Schwarz criterion (MSC):
    % msc(i) = logdp(i)/m - (log(ne) - 2.5) * (1 - 2.5*np(i)/(ne-np(i)));
    
    i      = i-1;                % go to next lower order
end

end



% y = Filter0(b, x)
%
% filters x with a fir filter so it has zero phase, i.e. shifts the
% filtered signal to the right half of the length of b.
%
% for now it zero pads the original signal
% later it might also do reflecton boundary conditions.
%
% be careful about the order of b!
% for even filters it is not exact  (change of Anton)
% - tired that even filterss dont' work


function y = Filter0In(b, x)

if size(x,1) == 1
	x = x(:);
end

% if mod(length(b),2)~=1
% 	error('filter order should be odd');
% end
if mod(length(b),2)~=1
    shift = length(b)/2;
else
    shift = (length(b)-1)/2;
end

[y0 z] = filter(b,1,single(x));

y = [y0(shift+1:end,:) ; z(1:shift,:)];

end


% computes the available memory in bytes
function HowMuch = FreeMemoryIn
if isunix
	[junk mem] = unix('vmstat |tail -1|awk ''{print $4} {print $6}''');
	HowMuch = sum(mem);
else
	HowMuch = 200;
	%200Mb for windows machin
	
end
end

function bool = inttoboolIn(ints,totalpoints)
% Takes a series of start-stop intervals (one row for each int, each row is
% [startpoint stoppoint]) and converts to a boolean with length
% (totalpoints) with zeros by default but 1s wherever points are inside
% "ints".  If totalpoints is not input then length is set by the last point
% in the last int.

warning off

if ~exist('totalpoints','var')
    totalpoints = 1;
end

bool = zeros(1,totalpoints);
for a = 1:size(ints,1);
    bool(round(ints(a,1)):round(ints(a,2))) = 1;
end
end


function FO = ViewAutoScoreThresholds(obj,event)
%Pull up window for user to look at where auto-scoring thresholds were and 
%also allow them to change them
FO = guidata(obj);
baseName = FO.baseName;
basePath = FO.basePath;

if ~isfield(FO,'AutoScore')
    FO.AutoScore = [];
end

% auto-load DetectionParameters, if they exist.  Save for later
SleepState = bz_LoadStates(basePath,'SleepState');

paramsAvailBool = 0;
detectnameAvailBool = 0;
if isempty(SleepState)
    answer = questdlg({'No SleepState.states.mat. Would you like to run SleepScoreMaster?', 'WARNING: will lose current states!!!'},...
        'AutoScore?');
    switch answer
        case 'Yes'
            %answer = questdlg('Use Loaded Channels?','Yes','No, Auto select channgels for scoring');
            SleepScoreMaster(basePath)
        case {'No','Cancel'}
            return
    end
elseif isfield(SleepState,'detectorinfo')
    if isfield(SleepState.detectorinfo,'detectorname')
        detectorname = SleepState.detectorinfo.detectorname;
        detectnameAvailBool = 1;
    end
    if isfield(SleepState.detectorinfo,'detectionparms')
        detectionparms = SleepState.detectorinfo.detectionparms;
        paramsAvailBool = 1;
    end
elseif isfield(SleepState,'detectorname') %old version of data
    detectorname = SleepState.detectorname;
    detectnameAvailBool = 1;
    if isfield(SleepState,'detectorparms')
        detectionparms = SleepState.detectorparms;
        paramsAvailBool = 1;
    end
end

if detectnameAvailBool
    switch detectorname
        case {'TheStateEditor'}
            answer = questdlg({'States manually detected. Would you like to run SleepScoreMaster?',...
                'WARNING: will lose current states!!!'},...
                'AutoScore?');
            switch answer
                case 'Yes'
                    SleepScoreMaster(basePath)
                case {'No','Cancel'}
                    return
            end
    end
else
    disp('No Detector name recorded or found')
end
        

if paramsAvailBool 
    %DL: this seems weird to me... 
    %why do we assume the detectionparms in SleepState are auto?
    FO.AutoScore.detectionparms =  detectionparms;
else
    FO.AutoScore.detectionparms = [];
end


% get histograms and thresholds that were used to determine states
HistAndThreshAlready_Bool = 0;
if isfield(FO,'AutoScore')
    if isfield(FO.AutoScore,'histsandthreshs')
        HistAndThreshAlready_Bool = 1;
    end
end
if HistAndThreshAlready_Bool
    histsandthreshs = FO.AutoScore.histsandthreshs;
else
    histsandthreshs = SSHistogramsAndThresholds_In(baseName,basePath);
    if isempty(histsandthreshs)
        disp('Exiting, no HistsAndThreshs')
        return
    end
end
if ~isfield(histsandthreshs,'stickySW'); histsandthreshs.stickySW = false; end
if ~isfield(histsandthreshs,'stickyTH'); histsandthreshs.stickyTH = false; end
if ~isfield(histsandthreshs,'stickyEMG');histsandthreshs.stickyEMG = false; end
FO.AutoScore.histsandthreshs = histsandthreshs;

% get histograms and thresholds of original detection
HistAndThreshOrigAlready_Bool = 0;
if isfield(SleepState,'detectorinfo')
    if isfield(SleepState.detectorinfo,'detectionparms')
        if isfield(SleepState.detectorinfo.detectionparms,'SleepScoreMetrics')
            if isfield(SleepState.detectorinfo.detectionparms.SleepScoreMetrics,'histsandthreshs_orig')
                HistAndThreshOrigAlready_Bool = 1;
            end
        end
    end
end
if HistAndThreshOrigAlready_Bool
    histsandthreshs_orig = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs_orig;
else
    histsandthreshs_orig = histsandthreshs;
end
FO.AutoScore.histsandthreshs_orig = histsandthreshs_orig;



% start figure
h = figure('position',[940 5 960 720]);
set(h, 'MenuBar', 'none');
set(h, 'ToolBar', 'none');

%top plot: Slow wave power
ax1 = subplot(3,2,2,'ButtonDownFcn',@ClickSetsLineXIn);hold on;
bar(histsandthreshs.swhistbins,histsandthreshs.swhist)
swline = plot(ax1,[histsandthreshs.swthresh histsandthreshs.swthresh],ylim(ax1),'r','tag','bw');
xlabel('SWS Band Power (NREM vs other)')
ylabel('Counts (sec)')
ResetToInitButton_sw = uicontrol('style', 'pushbutton', 'String', 'Init', 'Units', 'normalized', 'Position',  [0.92, 0.87, 0.08, 0.05]);
set(ResetToInitButton_sw,'Callback',@ResetToInitSw);
ResetToOrigButton_sw = uicontrol('style', 'pushbutton', 'String', 'Orig', 'Units', 'normalized', 'Position',  [0.92, 0.82, 0.08, 0.05]);
set(ResetToOrigButton_sw,'Callback',@ResetToOrigSw);
StickyThreshCheck_sw = uicontrol('style', 'checkbox', 'String', 'Sticky', 'Units', 'normalized', 'Position',  [0.92, 0.77, 0.08, 0.05],...
    'Value',histsandthreshs.stickySW);
%set(StickyThreshCheck_sw,'Callback',@SetStickySw);

title({'Click in plots to reset X value of thresholds',...
    'Setting thresholds to ''Sticky'' will reduce noise'})

%middle plot: EMG amplitude
ax2 = subplot(3,2,4,'ButtonDownFcn',@ClickSetsLineXIn);hold on;
bar(histsandthreshs.EMGhistbins,histsandthreshs.EMGhist)
EMGline = plot(ax2,[histsandthreshs.EMGthresh histsandthreshs.EMGthresh],ylim(ax2),'r','tag','bw');
xlabel('EMG (300-600Hz Correlation, active WAKE vs REM/inactive)')
ylabel('Counts (sec)')
ResetToInitButton_EMG = uicontrol('style', 'pushbutton', 'String', 'Init', 'Units', 'normalized', 'Position',  [0.92, 0.58, 0.08, 0.05]);
set(ResetToInitButton_EMG,'Callback',@ResetToInitEMG);
ResetToOrigButton_EMG = uicontrol('style', 'pushbutton', 'String', 'Orig', 'Units', 'normalized', 'Position',  [0.92, 0.53, 0.08, 0.05]);
set(ResetToOrigButton_EMG,'Callback',@ResetToOrigEMG);
StickyThreshCheck_EMG = uicontrol('style', 'checkbox', 'String', 'Sticky', 'Units', 'normalized', 'Position',  [0.92, 0.48, 0.08, 0.05],...
    'Value',histsandthreshs.stickyEMG);
%set(StickyThreshCheck_EMG,'Callback',@SetStickyEMG);

%bottom plot: Theta power
ax3 = subplot(3,2,6,'ButtonDownFcn',@ClickSetsLineXIn);hold on;
bar(histsandthreshs.THhistbins,histsandthreshs.THhist)
THline = plot(ax3,[histsandthreshs.THthresh histsandthreshs.THthresh],ylim(ax3),'r','tag','bw');
xlabel('Theta ratio (5-10Hz/2-15Hz, REM vs inactive WAKE)')
ylabel('Counts (sec)')
ResetToInitButton_TH = uicontrol('style', 'pushbutton', 'String', 'Init', 'Units', 'normalized', 'Position',  [0.92, 0.28, 0.08, 0.05]);
set(ResetToInitButton_TH,'Callback',@ResetToInitTH);
ResetToOrigButton_TH = uicontrol('style', 'pushbutton', 'String', 'Orig', 'Units', 'normalized', 'Position',  [0.92, 0.23, 0.08, 0.05]);
set(ResetToOrigButton_TH,'Callback',@ResetToOrigTH);
StickyThreshCheck_TH = uicontrol('style', 'checkbox', 'String', 'Sticky', 'Units', 'normalized', 'Position',  [0.92, 0.18, 0.08, 0.05],...
    'Value',histsandthreshs.stickyTH);
%set(StickyThreshCheck_TH,'Callback',@SetStickyTH);

%For 2d cluster plots
broadbandSlowWave = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.broadbandSlowWave;
thratio = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;
EMG = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.EMG;
plotstates = interp1(FO.to,FO.States,SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus,'nearest');
%right plot
ax4 = subplot(2,2,1);hold on;
    for ss = 1:5
        plot(ax4,broadbandSlowWave(plotstates==ss),EMG(plotstates==ss),'.','color',FO.colors.states{ss},'markersize',1)
    end
    swline2 = plot(histsandthreshs.swthresh*[1 1],ylim(ax4),'r','LineWidth',1);
    EMGline2 = plot(histsandthreshs.swthresh*[0 1],histsandthreshs.EMGthresh*[1 1],'r','LineWidth',1);
    xlabel('Broadband SW');ylabel('EMG')
    title('Isolate NREM')
    
    
ax5 = subplot(2,2,3);hold on;
    for ss = [1 2 4 5]
    plot(ax5,thratio(plotstates==ss),EMG(plotstates==ss),'.','color',FO.colors.states{ss},'markersize',1)
    end
    xlabel('Narrowband Theta');ylabel('EMG')
    thline2 = plot(histsandthreshs.THthresh*[1 1],histsandthreshs.EMGthresh*[0 1],'r','LineWidth',1);
    EMGline3 = plot([0 1],histsandthreshs.EMGthresh*[1 1],'r','LineWidth',1);
    title('Isolate REM from WAKE')

%RESCORE!
ReScoreButton = uicontrol('style', 'pushbutton', 'String', 'Re-Score', 'Units', 'normalized', 'Position',  [0.4, 0.01, 0.2, 0.05]);
set(ReScoreButton,'Callback',@ReClusterStates_In);


AutoClusterFig.fig = h;
AutoClusterFig.ax1 = ax1;
AutoClusterFig.ax2 = ax2;
AutoClusterFig.ax3 = ax3;
AutoClusterFig.ax4 = ax4;
AutoClusterFig.ax5 = ax5;
AutoClusterFig.swline = swline;
AutoClusterFig.EMGline = EMGline;
AutoClusterFig.THline = THline;
AutoClusterFig.swline2 = swline2;
AutoClusterFig.EMGline2 = EMGline2;
AutoClusterFig.thline2 = thline2;
AutoClusterFig.EMGline3 = EMGline3;
AutoClusterFig.stickySWbox = StickyThreshCheck_sw;
AutoClusterFig.stickyEMGbox = StickyThreshCheck_EMG;
AutoClusterFig.stickyTHbox = StickyThreshCheck_TH;
AutoClusterFig.histsandthreshs_init = histsandthreshs;%store first value

FO.AutoClusterFig = AutoClusterFig;
FO.AutoScore.histsandthreshs = histsandthreshs;

bool = 0;
if isfield(FO,'AutoScore')
    if isfield(FO.AutoScore,'histsandthreshs_orig')
        bool = 1;
    end
end
if ~bool
    FO.AutoScore.histsandthreshs_orig = histsandthreshs;
end

% guidata update is done after output
guidata(FO.fig,FO)

end


function ClickSetsLineXIn(obj,ev)
%For Autoscore fig (a push), lets a click set the location of a thresh line
cp = get(obj,'CurrentPoint');
newx = cp(1);
lo = findobj('parent',obj,'type','line','tag','bw');
set(lo,'XData',[newx newx])

end


function histsandthreshs = SSHistogramsAndThresholds_In(baseName,basePath)
% Get initial histograms and thresholds as calculated by SleepScoreMaster.m
SleepState = bz_LoadStates(basePath,'SleepState');

if isempty(SleepState) %If there is no saved SleepState already. 
    %Will also need something here if the name of the Detector isn't
    %'SleepScoreMaster'... sorry I didn't fix this -DL
    error('No SleepState.states.mat detected');
end

if isfield(SleepState,'detectorinfo')%handing different historical versions - for back compatibility
    dp = SleepState.detectorinfo.detectionparms;%current version
elseif isfield(SleepState,'detectorparams')
    dp = SleepState.detectorparams;%old version
end
    
%if already exists, just take from saved data
histsandthreshsOK = 0;
if isfield(dp,'histsandthreshs')
    histsandthreshs = dp.histsandthreshs;
    histsandthreshsOK = 1;
elseif isfield(dp,'SleepScoreMetrics')
    if isfield(dp.SleepScoreMetrics,'histsandthreshs')
        %this is the proper formatting, everything else in here is to deal
        %with legacy issues (DL 3/18/18)
        histsandthreshs = dp.SleepScoreMetrics.histsandthreshs;
    	histsandthreshsOK = 1;  
    end
end
%if not, try to re-make it
if ~histsandthreshsOK
    warning('We were unable to find histsandthreshs in your SleepState. Trying to recalculate...')

    loadgood = 1;
    if exist(fullfile(basePath,[baseName '.SleepScoreLFP.LFP.mat']),'file')
        s = load(fullfile(basePath,[baseName '.SleepScoreLFP.LFP.mat']));
    else
        loadgood = 0;%signify couldn't find preprocessed data
    end
    
    if exist(fullfile(basePath,[baseName '.EMGFromLFP.LFP.mat']),'file')
        e = load(fullfile(basePath,[baseName '.EMGFromLFP.LFP.mat']));
    else
        loadgood = 0;%signify couldn't find preprocessed data
    end
    
    
    if ~loadgood
        answer = questdlg({'No HistsAndThreshs.  Also no SleepscoreLFP.LFP.mat or .EMGFromLFP.LFP.mat found. Would you like to run SleepScoreMaster?', 'WARNING: Will lose current states!!!'},...
        'AutoScore?');
        switch answer
            case 'Yes'
                %answer = questdlg('Use Loaded Channels?','Yes','No, Auto select channgels for scoring');
                SleepScoreMaster(basePath)
                s = load(fullfile(basePath,[baseName '.SleepScoreLFP.LFP.mat']));
                e = load(fullfile(basePath,[baseName '.EMGFromLFP.LFP.mat']));
            case {'No','Cancel'}
                histsandthreshs = [];
                return
        end
    end

    
    [SleepScoreMetrics,StatePlotMaterials] = ClusterStates_GetMetrics(...
                                           basePath,s.SleepScoreLFP,e.EMGFromLFP,false);

    histsandthreshs = SleepScoreMetrics.histsandthreshs;

    SleepState.detectorinfo.detectionparms.SleepScoreMetrics = SleepScoreMetrics;
    SleepState.detectorinfo.StatePlotMaterials = StatePlotMaterials;
    save(fullfile(basePath,[baseName '.SleepState.states.mat']),'SleepState');
end
end

function [states] = ReClusterStates_In(obj,ev)
% Wrapper around functions ClusterStates_DetermineStates and
% StatesToEpisodes.  Note only the more raw (but not totally raw) 
% SleepState output from StatesToEpisodes is used... not the Episodes.
% One could consider excluding the refining step of StatesToEpidodes,
% but I think Dan Levenstein would not stand by that approach as
% appropriate and vetted
%

obj = findobj('tag','StateEditorMaster');
FO = guidata(obj(end));
baseName = FO.baseName;
basePath = FO.basePath;

% load detectionparameters if not loaded when user pressed "a" in ViewAutoScoreThresholds
if isfield(FO,'AutoScore')
    if isfield(FO.AutoScore,'detectionparms')
        dp = FO.AutoScore.detectionparms;
    end
end
if ~exist('dp','var')
    SleepState = bz_LoadStates(basePath,'SleepState');
    dp = SleepState.detectorinfo.detectionparms;
end

% grab user-entered thresholds from GUI, for final input to DetermineStates
swthresh = get(FO.AutoClusterFig.swline,'XData');
swthresh = swthresh(1,1);
EMGthresh = get(FO.AutoClusterFig.EMGline,'XData');
EMGthresh = EMGthresh(1,1);
THthresh = get(FO.AutoClusterFig.THline,'XData');
THthresh = THthresh(1,1);
FO.AutoScore.histsandthreshs.swthresh = swthresh;
FO.AutoScore.histsandthreshs.EMGthresh = EMGthresh;
FO.AutoScore.histsandthreshs.THthresh = THthresh;

FO.AutoScore.histsandthreshs.stickySW = FO.AutoClusterFig.stickySWbox.Value;
FO.AutoScore.histsandthreshs.stickyTH = FO.AutoClusterFig.stickyTHbox.Value;
FO.AutoScore.histsandthreshs.stickyEMG = FO.AutoClusterFig.stickyEMGbox.Value;


if ~isfield(dp,'MinTimeWindowParms')
    display('No MinTimeWindowParms found... using defaults')
    dp.MinTimeWindowParms = [];
end

% Execute scoring - USE SleepScore toolbox functions
[stateintervals,stateidx,~] = ClusterStates_DetermineStates(...
                                           dp.SleepScoreMetrics,dp.MinTimeWindowParms,FO.AutoScore.histsandthreshs);

plotstates = interp1(stateidx.timestamps,stateidx.states,dp.SleepScoreMetrics.t_clus,'nearest');
%Redraw 2d cluster plots
cla(FO.AutoClusterFig.ax4)
    for ss = 1:5
        plot(FO.AutoClusterFig.ax4,dp.SleepScoreMetrics.broadbandSlowWave(plotstates==ss),...
            dp.SleepScoreMetrics.EMG(plotstates==ss),'.','color',FO.colors.states{ss},'markersize',1)
    end
    plot(FO.AutoClusterFig.ax4,swthresh.*[1 1],ylim(FO.AutoClusterFig.ax4),'r','LineWidth',1);
    plot(FO.AutoClusterFig.ax4,[0 swthresh],EMGthresh.*[1 1],'r','LineWidth',1);
    
cla(FO.AutoClusterFig.ax5)
    for ss = [1 2 4 5]
        plot(FO.AutoClusterFig.ax5,dp.SleepScoreMetrics.thratio(plotstates==ss),...
            dp.SleepScoreMetrics.EMG(plotstates==ss),'.','color',FO.colors.states{ss},'markersize',1)
    end
    plot(FO.AutoClusterFig.ax5,THthresh.*[1 1],[0 EMGthresh],'r','LineWidth',1);
    plot(FO.AutoClusterFig.ax5,[0 1],EMGthresh.*[1 1],'r','LineWidth',1); 

%Interpolate to match fspec{1}.to
states = interp1(stateidx.timestamps,stateidx.states,FO.to,'nearest');
states(isnan(states)) = 0; %Bug fix DL - FO.to starts at time 0....    
    
FO.States = states;
guidata(FO.fig,FO);
modifyStates(1, states, 0);

end


function ResetToOrigSw(obj,ev)
obj = findobj('tag','StateEditorMaster');
FO = guidata(obj(end));
baseName = FO.baseName;
y = [0 max(FO.AutoScore.histsandthreshs_orig.swhist)];
x = [FO.AutoScore.histsandthreshs_orig.swthresh FO.AutoScore.histsandthreshs_orig.swthresh];
set(FO.AutoClusterFig.swline,'XData',x);
end

function ResetToOrigEMG(obj,ev)
obj = findobj('tag','StateEditorMaster');
FO = guidata(obj(end));
baseName = FO.baseName;
y = [0 max(FO.AutoScore.histsandthreshs_orig.EMGhist)];
x = [FO.AutoScore.histsandthreshs_orig.EMGthresh FO.AutoScore.histsandthreshs_orig.EMGthresh];
set(FO.AutoClusterFig.EMGline,'XData',x);
end

function ResetToOrigTH(obj,ev)
obj = findobj('tag','StateEditorMaster');
FO = guidata(obj(end));
baseName = FO.baseName;
y = [0 max(FO.AutoScore.histsandthreshs_orig.THhist)];
x = [FO.AutoScore.histsandthreshs_orig.THthresh FO.AutoScore.histsandthreshs_orig.THthresh];
set(FO.AutoClusterFig.THline,'XData',x);
end

function ResetToInitSw(obj,ev)
obj = findobj('tag','StateEditorMaster');
FO = guidata(obj(end));
baseName = FO.baseName;
y = [0 max(FO.AutoClusterFig.histsandthreshs_init.swhist)];
x = [FO.AutoClusterFig.histsandthreshs_init.swthresh FO.AutoClusterFig.histsandthreshs_init.swthresh];
set(FO.AutoClusterFig.swline,'XData',x);
end

function ResetToInitEMG(obj,ev)
obj = findobj('tag','StateEditorMaster');
FO = guidata(obj(end));
baseName = FO.baseName;
y = [0 max(FO.AutoClusterFig.histsandthreshs_init.EMGhist)];
x = [FO.AutoClusterFig.histsandthreshs_init.EMGthresh FO.AutoClusterFig.histsandthreshs_init.EMGthresh];
set(FO.AutoClusterFig.EMGline,'XData',x);
end

function ResetToInitTH(obj,ev)
obj = findobj('tag','StateEditorMaster');
FO = guidata(obj(end));
baseName = FO.baseName;
y = [0 max(FO.AutoClusterFig.histsandthreshs_init.THhist)];
x = [FO.AutoClusterFig.histsandthreshs_init.THthresh FO.AutoClusterFig.histsandthreshs_init.THthresh];
set(FO.AutoClusterFig.THline,'XData',x);
end

% function SetStickySw(obj,ev)
% obj = findobj('tag','StateEditorMaster');
% FO = guidata(obj(end));
% 
% set(FO.AutoClusterFig.swline,'XData',x);
% end

function [ INT ] = IDXtoINT_In( IDX ,numstates)
%IDXtoINT_In(IDX) Converts state indices to state on/offsets
%
%INPUT
%   IDX:    [t x 1] vector of state indices, where states are identified by
%           integers starting from 1. Times with IDX 0 will not be counted
%           in any interval INT
%   numstates (optional)  number of interval types (for use
%
%OUTPUT
%   INT:    {nstates} cell array of intervals - start and end times
%
%DLevenstein 2015-16
%%
if ~exist('numstates','var')
    numstates = max(IDX);
end
states = 1:numstates;
if isrow(IDX)
    IDX = IDX';
end
IDX = [0; IDX; 0];
for ss = 1:numstates
    statetimes = IDX==states(ss);
    INT{ss} = [find(diff(statetimes)==1) find(diff(statetimes)==-1)-1];
end

end


function newvals = ResampleTolerant_IN(vals,length1,length2)
% Wrapper around the resample function that allows it to work even if
% the product of the lengths are long enough to overwhelm the resample.m 
% limit of 2^31
% Works by finding a rational number approximation of the requested length
% ratio... to within a particular tolerance.
% INPUTS
% vals = vector of values to be resampled
% length1 = desired length
% length2 = initial length (often equals length(vals))
%
% Dan Levenstein code made into a function by Brendon Watson
% August 2016


if length1*length2 < 2^31 % if no need to change factors don't
    newvals = resample(vals,length1,length2);
else
    newvals = [1 1];
    
    resamplefact = length1/length2;
    tol = 0.0001;
    while length(newvals(:,1)) ~= length1
        [P,Q] = rat(resamplefact,tol);
        if P==0
            tol = tol/10;
            continue
        end
        if P*Q >=2^20 || tol<1e-300  %Avoid crashing resample...
            vals([1,end],:) = [];
            length2 = length2-2;
            resamplefact = length1/length2;
            tol = 0.0001;
            continue
        end
        newvals = resample(vals,P,Q);
        tol = tol/10;
    end
end

end


function [Ypk,Xpk,Wpk,Ppk] = findpeaks_In(Yin,varargin)
%FINDPEAKS Find local peaks in data
%   PKS = FINDPEAKS(Y) finds local peaks in the data vector Y. A local peak
%   is defined as a data sample which is either larger than the two
%   neighboring samples or is equal to Inf.
%
%   [PKS,LOCS]= FINDPEAKS(Y) also returns the indices LOCS at which the
%   peaks occur.
%
%   [PKS,LOCS] = FINDPEAKS(Y,X) specifies X as the location vector of data
%   vector Y. X must be a strictly increasing vector of the same length as
%   Y. LOCS returns the corresponding value of X for each peak detected.
%   If X is omitted, then X will correspond to the indices of Y.
%
%   [PKS,LOCS] = FINDPEAKS(Y,Fs) specifies the sample rate, Fs, as a
%   positive scalar, where the first sample instant of Y corresponds to a
%   time of zero.
%
%   [...] = FINDPEAKS(...,'MinPeakHeight',MPH) finds only those peaks that
%   are greater than the minimum peak height, MPH. MPH is a real valued
%   scalar. The default value of MPH is -Inf.
%
%   [...] = FINDPEAKS(...,'MinPeakProminence',MPP) finds peaks guaranteed
%   to have a vertical drop of more than MPP from the peak on both sides
%   without encountering either the end of the signal or a larger
%   intervening peak. The default value of MPP is zero.
%
%   [...] = FINDPEAKS(...,'Threshold',TH) finds peaks that are at least
%   greater than both adjacent samples by the threshold, TH. TH is real
%   valued scalar greater than or equal to zero. The default value of TH is
%   zero.
%
%   FINDPEAKS(...,'WidthReference',WR) estimates the width of the peak as
%   the distance between the points where the signal intercepts a
%   horizontal reference line. The points are found by linear
%   interpolation. The height of the line is selected using the criterion
%   specified in WR:
% 
%    'halfprom' - the reference line is positioned beneath the peak at a
%       vertical distance equal to half the peak prominence.
% 
%    'halfheight' - the reference line is positioned at one-half the peak 
%       height. The line is truncated if any of its intercept points lie
%       beyond the borders of the peaks selected by the 'MinPeakHeight',
%       'MinPeakProminence' and 'Threshold' parameters. The border between
%       peaks is defined by the horizontal position of the lowest valley
%       between them. Peaks with heights less than zero are discarded.
% 
%    The default value of WR is 'halfprom'.
%
%   [...] = FINDPEAKS(...,'MinPeakWidth',MINW) finds peaks whose width is
%   at least MINW. The default value of MINW is zero.
%
%   [...] = FINDPEAKS(...,'MaxPeakWidth',MAXW) finds peaks whose width is
%   at most MAXW. The default value of MAXW is Inf.
%
%   [...] = FINDPEAKS(...,'MinPeakDistance',MPD) finds peaks separated by
%   more than the minimum peak distance, MPD. This parameter may be
%   specified to ignore smaller peaks that may occur in close proximity to
%   a large local peak. For example, if a large local peak occurs at LOC,
%   then all smaller peaks in the range [N-MPD, N+MPD] are ignored. If not
%   specified, MPD is assigned a value of zero.
%
%   [...] = FINDPEAKS(...,'SortStr',DIR) specifies the direction of sorting
%   of peaks. DIR can take values of 'ascend', 'descend' or 'none'. If not
%   specified, DIR takes the value of 'none' and the peaks are returned in
%   the order of their occurrence.
%
%   [...] = FINDPEAKS(...,'NPeaks',NP) specifies the maximum number of peaks
%   to be found. NP is an integer greater than zero. If not specified, all
%   peaks are returned. Use this parameter in conjunction with setting the
%   sort direction to 'descend' to return the NP largest peaks. (see
%   'SortStr')
%
%   [PKS,LOCS,W] = FINDPEAKS(...) returns the width, W, of each peak by
%   linear interpolation of the left- and right- intercept points to the
%   reference defined by 'WidthReference'.
%
%   [PKS,LOCS,W,P] = FINDPEAKS(...) returns the prominence, P, of each
%   peak.
%
%   FINDPEAKS(...) without output arguments plots the signal and the peak
%   values it finds
%
%   FINDPEAKS(...,'Annotate',PLOTSTYLE) will annotate a plot of the
%   signal with PLOTSTYLE. If PLOTSTYLE is 'peaks' the peaks will be
%   plotted. If PLOTSTYLE is 'extents' the signal, peak values, widths,
%   prominences of each peak will be annotated. 'Annotate' will be ignored
%   if called with output arguments. The default value of PLOTSTYLE is
%   'peaks'.
%
%   % Example 1:
%   %   Plot the Zurich numbers of sunspot activity from years 1700-1987
%   %   and identify all local maxima at least six years apart
%   load sunspot.dat
%   findpeaks(sunspot(:,2),sunspot(:,1),'MinPeakDistance',6)
%   xlabel('Year');
%   ylabel('Zurich number');
%
%   % Example 2: 
%   %   Plot peak values of an audio signal that drop at least 1V on either
%   %   side without encountering values larger than the peak.
%   load mtlb
%   findpeaks(mtlb,Fs,'MinPeakProminence',1)
%
%   % Example 3:
%   %   Plot all peaks of a chirp signal whose widths are between .5 and 1 
%   %   milliseconds.
%   Fs = 44.1e3; N = 1000;
%   x = sin(2*pi*(1:N)/N + (10*(1:N)/N).^2);
%   findpeaks(x,Fs,'MinPeakWidth',.5e-3,'MaxPeakWidth',1e-3, ...
%             'Annotate','extents')

%   Copyright 2007-2014 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

cond = nargin >= 1;
if ~cond
    coder.internal.assert(cond,'MATLAB:narginchk:notEnoughInputs');
end

cond = nargin <= 22;
if ~cond
    coder.internal.assert(cond,'MATLAB:narginchk:tooManyInputs');
end

% extract the parameters from the input argument list
[y,yIsRow,x,xIsRow,minH,minP,minW,maxW,minD,minT,maxN,sortDir,annotate,refW] ...
  = parse_inputs(Yin,varargin{:});

% find indices of all finite and infinite peaks and the inflection points
[iFinite,iInfite,iInflect] = getAllPeaks(y);

% keep only the indices of finite peaks that meet the required 
% minimum height and threshold
iPk = removePeaksBelowMinPeakHeight(y,iFinite,minH,refW);
iPk = removePeaksBelowThreshold(y,iPk,minT);

% indicate if we need to compute the extent of a peak
needWidth = minW>0 || maxW<inf || minP>0 || nargout>2 || strcmp(annotate,'extents');

if needWidth
  % obtain the indices of each peak (iPk), the prominence base (bPk), and
  % the x- and y- coordinates of the peak base (bxPk, byPk) and the width
  % (wxPk)
  [iPk,bPk,bxPk,byPk,wxPk] = findExtents(y,x,iPk,iFinite,iInfite,iInflect,minP,minW,maxW,refW);
else
  % combine finite and infinite peaks into one list
  [iPk,bPk,bxPk,byPk,wxPk] = combinePeaks(iPk,iInfite);
end

% find the indices of the largest peaks within the specified distance
idx = findPeaksSeparatedByMoreThanMinPeakDistance(y,x,iPk,minD);

% re-order and bound the number of peaks based upon the index vector
idx = orderPeaks(y,iPk,idx,sortDir);
idx = keepAtMostNpPeaks(idx,maxN);

% use the index vector to fetch the correct peaks.
iPk = iPk(idx);
if needWidth
  [bPk, bxPk, byPk, wxPk] = fetchPeakExtents(idx,bPk,bxPk,byPk,wxPk);
end

if nargout > 0
  % assign output variables
  if needWidth
    [Ypk,Xpk,Wpk,Ppk] = assignFullOutputs(y,x,iPk,wxPk,bPk,yIsRow,xIsRow);
  else
    [Ypk,Xpk] = assignOutputs(y,x,iPk,yIsRow,xIsRow);
  end    
else
  % no output arguments specified. plot and optionally annotate
  hAxes = plotSignalWithPeaks(x,y,iPk);
  if strcmp(annotate,'extents')
    plotExtents(hAxes,x,y,iPk,bPk,bxPk,byPk,wxPk,refW);
  end
  
  scalePlot(hAxes);
end

end
%--------------------------------------------------------------------------
function [y,yIsRow,x,xIsRow,Ph,Pp,Wmin,Wmax,Pd,Th,NpOut,Str,Ann,Ref] = parse_inputs(Yin,varargin)

% Validate input signal
validateattributes(Yin,{'numeric'},{'nonempty','real','vector'},...
    'findpeaks','Y');
yIsRow = isrow(Yin);
y = Yin(:);

% copy over orientation of y to x.
xIsRow = yIsRow;

% indicate if the user specified an Fs or X
hasX = ~isempty(varargin) && isnumeric(varargin{1});

if hasX
  startArg = 2;
  if isscalar(varargin{1})
    % Fs
    Fs = varargin{1};
    validateattributes(Fs,{'double'},{'real','finite','positive'},'findpeaks','Fs');
    x = (0:numel(y)-1).'/Fs;
  else
    % X
    Xin = varargin{1};
    validateattributes(Xin,{'double'},{'real','finite','vector','increasing'},'findpeaks','X');
    if numel(Xin) ~= numel(Yin)
      if coder.target('MATLAB')
        throwAsCaller(MException(message('signal:findpeaks:mismatchYX')));
      else
        coder.internal.errorIf(true,'signal:findpeaks:mismatchYX');
      end
    end
    xIsRow = isrow(Xin);
    x = Xin(:);
  end
else
  startArg = 1;
  % unspecified, use index vector
  x = (1:numel(y)).';
end

if coder.target('MATLAB')
    try %#ok<EMTC>
        % Check the input data type. Single precision is not supported.
        chkinputdatatype_In(y);
        chkinputdatatype_In(x);
    catch ME
        throwAsCaller(ME);
    end
else
    chkinputdatatype_In(y);
    chkinputdatatype_In(x);
end

M = numel(y);
cond = (M < 3);
if cond
    coder.internal.errorIf(cond,'signal:findpeaks:emptyDataSet');
end

%#function dspopts.findpeaks
defaultMinPeakHeight = -inf;
defaultMinPeakProminence = 0;
defaultMinPeakWidth = 0;
defaultMaxPeakWidth = Inf;
defaultMinPeakDistance = 0;
defaultThreshold = 0;
defaultNPeaks = [];
defaultSortStr = 'none';
defaultAnnotate = 'peaks';
defaultWidthReference = 'halfprom';

if coder.target('MATLAB')
    p = inputParser;
    addParameter(p,'MinPeakHeight',defaultMinPeakHeight);
    addParameter(p,'MinPeakProminence',defaultMinPeakProminence);
    addParameter(p,'MinPeakWidth',defaultMinPeakWidth);
    addParameter(p,'MaxPeakWidth',defaultMaxPeakWidth);
    addParameter(p,'MinPeakDistance',defaultMinPeakDistance);
    addParameter(p,'Threshold',defaultThreshold);
    addParameter(p,'NPeaks',defaultNPeaks);
    addParameter(p,'SortStr',defaultSortStr);
    addParameter(p,'Annotate',defaultAnnotate);
    addParameter(p,'WidthReference',defaultWidthReference);
    parse(p,varargin{startArg:end});
    Ph = p.Results.MinPeakHeight;
    Pp = p.Results.MinPeakProminence;
    Wmin = p.Results.MinPeakWidth;
    Wmax = p.Results.MaxPeakWidth;
    Pd = p.Results.MinPeakDistance;
    Th = p.Results.Threshold;
    Np = p.Results.NPeaks;
    Str = p.Results.SortStr;
    Ann = p.Results.Annotate;
    Ref = p.Results.WidthReference;
else
    parms = struct('MinPeakHeight',uint32(0), ...
                'MinPeakProminence',uint32(0), ...
                'MinPeakWidth',uint32(0), ...
                'MaxPeakWidth',uint32(0), ...
                'MinPeakDistance',uint32(0), ...
                'Threshold',uint32(0), ...
                'NPeaks',uint32(0), ...
                'SortStr',uint32(0), ...
                'Annotate',uint32(0), ...
                'WidthReference',uint32(0));
    pstruct = eml_parse_parameter_inputs(parms,[],varargin{startArg:end});
    Ph = eml_get_parameter_value(pstruct.MinPeakHeight,defaultMinPeakHeight,varargin{startArg:end});
    Pp = eml_get_parameter_value(pstruct.MinPeakProminence,defaultMinPeakProminence,varargin{startArg:end});
    Wmin = eml_get_parameter_value(pstruct.MinPeakWidth,defaultMinPeakWidth,varargin{startArg:end});
    Wmax = eml_get_parameter_value(pstruct.MaxPeakWidth,defaultMaxPeakWidth,varargin{startArg:end});
    Pd = eml_get_parameter_value(pstruct.MinPeakDistance,defaultMinPeakDistance,varargin{startArg:end});
    Th = eml_get_parameter_value(pstruct.Threshold,defaultThreshold,varargin{startArg:end});
    Np = eml_get_parameter_value(pstruct.NPeaks,defaultNPeaks,varargin{startArg:end});
    Str = eml_get_parameter_value(pstruct.SortStr,defaultSortStr,varargin{startArg:end});
    Ann = eml_get_parameter_value(pstruct.Annotate,defaultAnnotate,varargin{startArg:end});
    Ref = eml_get_parameter_value(pstruct.WidthReference,defaultWidthReference,varargin{startArg:end});
end

% limit the number of peaks to the number of input samples
if isempty(Np)
    NpOut = M;
else
    NpOut = Np;
end

% ignore peaks below zero when using halfheight width reference
if strcmp(Ref,'halfheight')
  Ph = max(Ph,0);
end

validateattributes(Ph,{'numeric'},{'real','scalar','nonempty'},'findpeaks','MinPeakHeight');
validateattributes(Pd,{'numeric'},{'real','scalar','nonempty','nonnegative','<',x(M)-x(1)},'findpeaks','MinPeakDistance');
validateattributes(Pp,{'numeric'},{'real','scalar','nonempty','nonnegative'},'findpeaks','MinPeakProminence');
validateattributes(Wmin,{'numeric'},{'real','scalar','finite','nonempty','nonnegative'},'findpeaks','MinPeakWidth');
validateattributes(Wmax,{'numeric'},{'real','scalar','nonnan','nonempty','nonnegative'},'findpeaks','MaxPeakWidth');
validateattributes(Pd,{'numeric'},{'real','scalar','nonempty','nonnegative'},'findpeaks','MinPeakDistance');
validateattributes(Th,{'numeric'},{'real','scalar','nonempty','nonnegative'},'findpeaks','Threshold');
validateattributes(NpOut,{'numeric'},{'real','scalar','nonempty','integer','positive'},'findpeaks','NPeaks');
Str = validatestring(Str,{'ascend','none','descend'},'findpeaks','SortStr');
Ann = validatestring(Ann,{'peaks','extents'},'findpeaks','SortStr');
Ref = validatestring(Ref,{'halfprom','halfheight'},'findpeaks','WidthReference');

end
%--------------------------------------------------------------------------
function [iPk,iInf,iInflect] = getAllPeaks(y)
% fetch indices all infinite peaks
iInf = find(isinf(y) & y>0);

% temporarily remove all +Inf values
yTemp = y;
yTemp(iInf) = NaN;

% determine the peaks and inflection points of the signal
[iPk,iInflect] = findLocalMaxima(yTemp);


end
%--------------------------------------------------------------------------
function [iPk, iInflect] = findLocalMaxima(yTemp)
% bookend Y by NaN and make index vector
yTemp = [NaN; yTemp; NaN];
iTemp = (1:numel(yTemp)).';

% keep only the first of any adjacent pairs of equal values (including NaN).
yFinite = ~isnan(yTemp);
iNeq = [1; 1 + find((yTemp(1:end-1) ~= yTemp(2:end)) & ...
                    (yFinite(1:end-1) | yFinite(2:end)))];
iTemp = iTemp(iNeq);

% take the sign of the first sample derivative
s = sign(diff(yTemp(iTemp)));

% find local maxima
iMax = 1 + find(diff(s)<0);

% find all transitions from rising to falling or to NaN
iAny = 1 + find(s(1:end-1)~=s(2:end));

% index into the original index vector without the NaN bookend.
iInflect = iTemp(iAny)-1;
iPk = iTemp(iMax)-1;

end
%--------------------------------------------------------------------------
function iPk = removePeaksBelowMinPeakHeight(Y,iPk,Ph,widthRef)
if ~isempty(iPk) 
  iPk = iPk(Y(iPk) > Ph);
  if isempty(iPk) && ~strcmp(widthRef,'halfheight')
    if coder.target('MATLAB')
        warning(message('signal:findpeaks:largeMinPeakHeight', 'MinPeakHeight', 'MinPeakHeight'));
    end
  end
end

end
%--------------------------------------------------------------------------
function iPk = removePeaksBelowThreshold(Y,iPk,Th)

base = max(Y(iPk-1),Y(iPk+1));
iPk = iPk(Y(iPk)-base >= Th);

end
%--------------------------------------------------------------------------
function [iPk,bPk,bxPk,byPk,wxPk] = findExtents(y,x,iPk,iFin,iInf,iInflect,minP,minW,maxW,refW)
% temporarily filter out +Inf from the input
yFinite = y;
yFinite(iInf) = NaN;

% get the base and left and right indices of each prominence base
[bPk,iLB,iRB] = getPeakBase(yFinite,iPk,iFin,iInflect);

% keep only those indices with at least the specified prominence
[iPk,bPk,iLB,iRB] = removePeaksBelowMinPeakProminence(yFinite,iPk,bPk,iLB,iRB,minP);

% get the x-coordinates of the half-height width borders of each peak
[wxPk,iLBh,iRBh] = getPeakWidth(yFinite,x,iPk,bPk,iLB,iRB,refW);

% merge finite and infinite peaks together into one list
[iPk,bPk,bxPk,byPk,wxPk] = combineFullPeaks(y,x,iPk,bPk,iLBh,iRBh,wxPk,iInf);

% keep only those in the range minW < w < maxW
[iPk,bPk,bxPk,byPk,wxPk] = removePeaksOutsideWidth(iPk,bPk,bxPk,byPk,wxPk,minW,maxW);

end

%--------------------------------------------------------------------------
function [peakBase,iLeftSaddle,iRightSaddle] = getPeakBase(yTemp,iPk,iFin,iInflect)
% determine the indices that border each finite peak
[iLeftBase, iLeftSaddle] = getLeftBase(yTemp,iPk,iFin,iInflect);
[iRightBase, iRightSaddle] = getLeftBase(yTemp,flipud(iPk),flipud(iFin),flipud(iInflect));
iRightBase = flipud(iRightBase);
iRightSaddle = flipud(iRightSaddle);
peakBase = max(yTemp(iLeftBase),yTemp(iRightBase));

end
%--------------------------------------------------------------------------
function [iBase, iSaddle] = getLeftBase(yTemp,iPeak,iFinite,iInflect)
% pre-initialize output base and saddle indices
iBase = zeros(size(iPeak));
iSaddle = zeros(size(iPeak));

% table stores the most recently encountered peaks in order of height
peak = zeros(size(iFinite));
valley = zeros(size(iFinite));
iValley = zeros(size(iFinite));

n = 0;
i = 1;
j = 1;
k = 1;

% pre-initialize v for code generation
v = NaN; 
iv = 1;

while k<=numel(iPeak)
  % walk through the inflections until you reach a peak
  while iInflect(i) ~= iFinite(j) 
    v = yTemp(iInflect(i));
    iv = iInflect(i);
    if isnan(v)
      % border seen, start over.
      n = 0;
    else
      % ignore previously stored peaks with a valley larger than this one
      while n>0 && valley(n)>v
        n = n - 1;
      end
    end
    i = i + 1;
  end
  % get the peak
  p = yTemp(iInflect(i));
  
  % keep the smallest valley of all smaller peaks
  while n>0 && peak(n) < p
    if valley(n) < v
      v = valley(n);
      iv = iValley(n);
    end
    n = n - 1;
  end

  % record "saddle" valleys in between equal-height peaks
  isv = iv;
  
  % keep seeking smaller valleys until you reach a larger peak
  while n>0 && peak(n) <= p
    if valley(n) < v
      v = valley(n);
      iv = iValley(n);
    end
    n = n - 1;      
  end
  
  % record the new peak and save the index of the valley into the base
  % and saddle
  n = n + 1;
  peak(n) = p;
  valley(n) = v;
  iValley(n) = iv;

  if iInflect(i) == iPeak(k)
    iBase(k) = iv;
    iSaddle(k) = isv;
    k = k + 1;
  end
  
  i = i + 1;
  j = j + 1;
end

end
%--------------------------------------------------------------------------
function [iPk,pbPk,iLB,iRB] = removePeaksBelowMinPeakProminence(y,iPk,pbPk,iLB,iRB,minP)
% compute the prominence of each peak
Ppk = y(iPk)-pbPk;

% keep those that are above the specified prominence
idx = find(Ppk >= minP);
iPk = iPk(idx);
pbPk = pbPk(idx);
iLB = iLB(idx);
iRB = iRB(idx);

end
%--------------------------------------------------------------------------
function [wxPk,iLBh,iRBh] = getPeakWidth(y,x,iPk,pbPk,iLB,iRB,wRef)
if isempty(iPk)
  % no peaks.  define empty containers
  base = zeros(size(iPk));
  iLBh = zeros(size(iPk));
  iRBh = zeros(size(iPk));  
elseif strcmp(wRef,'halfheight')
  % set the baseline to zero
  base = zeros(size(iPk));

  % border the width by no more than the lowest valley between this peak
  % and the next peak
  iLBh = [iLB(1); max(iLB(2:end),iRB(1:end-1))];
  iRBh = [min(iRB(1:end-1),iLB(2:end)); iRB(end)];
  iGuard = iLBh > iPk;
  iLBh(iGuard) = iLB(iGuard);
  iGuard = iRBh < iPk;
  iRBh(iGuard) = iRB(iGuard);
else
  % use the prominence base
  base = pbPk;
  
  % border the width by the saddle of the peak
  iLBh = iLB;
  iRBh = iRB;
end

% get the width boundaries of each peak
wxPk = getHalfMaxBounds(y, x, iPk, base, iLBh, iRBh);

end
%--------------------------------------------------------------------------
function bounds = getHalfMaxBounds(y, x, iPk, base, iLB, iRB)
bounds = zeros(numel(iPk),2);

% interpolate both the left and right bounds clamping at borders
for i=1:numel(iPk)
  
  % compute the desired reference level at half-height or half-prominence
  refHeight = (y(iPk(i))+base(i))/2;
  
  % compute the index of the left-intercept at half max
  iLeft = findLeftIntercept(y, iPk(i), iLB(i), refHeight);
  if iLeft < iLB(i)
    xLeft = x(iLB(i));
  else
    xLeft = linterp(x(iLeft),x(iLeft+1),y(iLeft),y(iLeft+1),y(iPk(i)),base(i));
  end
  
  % compute the index of the right-intercept
  iRight = findRightIntercept(y, iPk(i), iRB(i), refHeight);
  if iRight > iRB(i)
    xRight = x(iRB(i));
  else
    xRight = linterp(x(iRight), x(iRight-1), y(iRight), y(iRight-1), y(iPk(i)),base(i));
  end

  % store result
  bounds(i,:) = [xLeft xRight];
end

end
%--------------------------------------------------------------------------
function idx = findLeftIntercept(y, idx, borderIdx, refHeight)
% decrement index until you pass under the reference height or pass the
% index of the left border, whichever comes first
while idx>=borderIdx && y(idx) > refHeight
  idx = idx - 1;
end

end
%--------------------------------------------------------------------------
function idx = findRightIntercept(y, idx, borderIdx, refHeight)
% increment index until you pass under the reference height or pass the
% index of the right border, whichever comes first
while idx<=borderIdx && y(idx) > refHeight
  idx = idx + 1;
end

end
%--------------------------------------------------------------------------
function xc = linterp(xa,xb,ya,yb,yc,bc)
% interpolate between points (xa,ya) and (xb,yb) to find (xc, 0.5*(yc-yc)).
xc = xa + (xb-xa) .* (0.5*(yc+bc)-ya) ./ (yb-ya);

% invoke L'Hospital's rule when -Inf is encountered. 
if isnan(xc)
  % yc and yb are guaranteed to be finite. 
  if isinf(bc)
    % both ya and bc are -Inf.
    xc = 0.5*(xa+xb);
  else
    % only ya is -Inf.
    xc = xb;
  end
end

end
%--------------------------------------------------------------------------
function [iPk,bPk,bxPk,byPk,wxPk] = removePeaksOutsideWidth(iPk,bPk,bxPk,byPk,wxPk,minW,maxW)

if isempty(iPk) || minW==0 && maxW == inf
  return
end

% compute the width of each peak and extract the matching indices
w = diff(wxPk,1,2);
idx = find(minW <= w & w <= maxW);

% fetch the surviving peaks
iPk = iPk(idx);
bPk = bPk(idx);
bxPk = bxPk(idx,:);
byPk = byPk(idx,:);
wxPk = wxPk(idx,:);

end
%--------------------------------------------------------------------------
function [iPkOut,bPk,bxPk,byPk,wxPk] = combinePeaks(iPk,iInf)
iPkOut = union(iPk,iInf);
bPk = zeros(0,1);
bxPk = zeros(0,2);
byPk = zeros(0,2);
wxPk = zeros(0,2);

end
%--------------------------------------------------------------------------
function [iPkOut,bPkOut,bxPkOut,byPkOut,wxPkOut] = combineFullPeaks(y,x,iPk,bPk,iLBw,iRBw,wPk,iInf)
iPkOut = union(iPk, iInf);

% create map of new indices to old indices
[~, iFinite] = intersect(iPkOut,iPk);
[~, iInfinite] = intersect(iPkOut,iInf);

% prevent row concatenation when iPk and iInf both have less than one
% element
iPkOut = iPkOut(:);

% compute prominence base
bPkOut = zeros(size(iPkOut));
bPkOut(iFinite) = bPk;
bPkOut(iInfinite) = 0;

% compute indices of left and right infinite borders
iInfL = max(1,iInf-1);
iInfR = min(iInf+1,numel(x));

% copy out x- values of the left and right prominence base
% set each base border of an infinite peaks halfway between itself and
% the next adjacent sample
bxPkOut = zeros(size(iPkOut,1),2);
bxPkOut(iFinite,1) = x(iLBw);
bxPkOut(iFinite,2) = x(iRBw);
bxPkOut(iInfinite,1) = 0.5*(x(iInf)+x(iInfL));
bxPkOut(iInfinite,2) = 0.5*(x(iInf)+x(iInfR));

% copy out y- values of the left and right prominence base
byPkOut = zeros(size(iPkOut,1),2);
byPkOut(iFinite,1) = y(iLBw);
byPkOut(iFinite,2) = y(iRBw);
byPkOut(iInfinite,1) = y(iInfL);
byPkOut(iInfinite,2) = y(iInfR);

% copy out x- values of the width borders
% set each width borders of an infinite peaks halfway between itself and
% the next adjacent sample
wxPkOut = zeros(size(iPkOut,1),2);
wxPkOut(iFinite,:) = wPk;
wxPkOut(iInfinite,1) = 0.5*(x(iInf)+x(iInfL));
wxPkOut(iInfinite,2) = 0.5*(x(iInf)+x(iInfR));

end
%--------------------------------------------------------------------------
function idx = findPeaksSeparatedByMoreThanMinPeakDistance(y,x,iPk,Pd)
% Start with the larger peaks to make sure we don't accidentally keep a
% small peak and remove a large peak in its neighborhood. 

if isempty(iPk) || Pd==0
  idx = (1:numel(iPk)).';
  return
end

% copy peak values and locations to a temporary place
pks = y(iPk);
locs = x(iPk);

% Order peaks from large to small
[~, sortIdx] = sort(pks,'descend');
locs_temp = locs(sortIdx);

idelete = ones(size(locs_temp))<0;
for i = 1:length(locs_temp)
  if ~idelete(i)
    % If the peak is not in the neighborhood of a larger peak, find
    % secondary peaks to eliminate.
    idelete = idelete | (locs_temp>=locs_temp(i)-Pd)&(locs_temp<=locs_temp(i)+Pd); 
    idelete(i) = 0; % Keep current peak
  end
end

% report back indices in consecutive order
idx = sort(sortIdx(~idelete));


end

%--------------------------------------------------------------------------
function idx = orderPeaks(Y,iPk,idx,Str)

if isempty(idx) || strcmp(Str,'none')
  return
end

if strcmp(Str,'ascend')
  [~,s]  = sort(Y(iPk(idx)),'ascend');
else
  [~,s]  = sort(Y(iPk(idx)),'descend');
end

idx = idx(s);


end
%--------------------------------------------------------------------------
function idx = keepAtMostNpPeaks(idx,Np)

if length(idx)>Np
  idx = idx(1:Np);
end

end
%--------------------------------------------------------------------------
function [bPk,bxPk,byPk,wxPk] = fetchPeakExtents(idx,bPk,bxPk,byPk,wxPk)
bPk = bPk(idx);
bxPk = bxPk(idx,:);
byPk = byPk(idx,:);
wxPk = wxPk(idx,:);

end
%--------------------------------------------------------------------------
function [YpkOut,XpkOut] = assignOutputs(y,x,iPk,yIsRow,xIsRow)

% fetch the coordinates of the peak
Ypk = y(iPk);
Xpk = x(iPk);

% preserve orientation of Y
if yIsRow
  YpkOut = Ypk.';
else
  YpkOut = Ypk;
end

% preserve orientation of X
if xIsRow
  XpkOut = Xpk.';
else
  XpkOut = Xpk;
end

end
%--------------------------------------------------------------------------
function [YpkOut,XpkOut,WpkOut,PpkOut] = assignFullOutputs(y,x,iPk,wxPk,bPk,yIsRow,xIsRow)

% fetch the coordinates of the peak
Ypk = y(iPk);
Xpk = x(iPk);

% compute the width and prominence
Wpk = diff(wxPk,1,2);
Ppk = Ypk-bPk;

% preserve orientation of Y (and P)
if yIsRow
  YpkOut = Ypk.';
  PpkOut = Ppk.';
else
  YpkOut = Ypk;
  PpkOut = Ppk;  
end

% preserve orientation of X (and W)
if xIsRow
  XpkOut = Xpk.';
  WpkOut = Wpk.';
else
  XpkOut = Xpk;
  WpkOut = Wpk;  
end

end
%--------------------------------------------------------------------------
function hAxes = plotSignalWithPeaks(x,y,iPk)

% plot signal
hLine = plot(x,y,'Tag','Signal');
hAxes = ancestor(hLine,'Axes');
% turn on grid
grid on;

% use the color of the line
color = get(hLine,'Color');
hLine = line(x(iPk),y(iPk),'Parent',hAxes, ...
     'Marker','o','LineStyle','none','Color',color,'tag','Peak');

% if using MATLAB use offset inverted triangular marker
if coder.target('MATLAB')
  plotpkmarkers(hLine,y(iPk));
end

end

%--------------------------------------------------------------------------
function plotExtents(hAxes,x,y,iPk,bPk,bxPk,byPk,wxPk,refW)

% compute level of half-maximum (height or prominence)
if strcmp(refW,'halfheight')
  hm = 0.5*y(iPk);
else
  hm = 0.5*(y(iPk)+bPk);
end

% get the default color order
colors = get(0,'DefaultAxesColorOrder');

% plot boundaries between adjacent peaks when using half-height
if strcmp(refW,'halfheight')
  % plot height
  plotLines(hAxes,'Height',x(iPk),y(iPk),x(iPk),zeros(numel(iPk),1),colors(2,:));  

  % plot width
  plotLines(hAxes,'HalfHeightWidth',wxPk(:,1),hm,wxPk(:,2),hm,colors(3,:));
      
  % plot peak borders
  idx = find(byPk(:,1)>0);
  plotLines(hAxes,'Border',bxPk(idx,1),zeros(numel(idx),1),bxPk(idx,1),byPk(idx,1),colors(4,:));
  idx = find(byPk(:,2)>0);
  plotLines(hAxes,'Border',bxPk(idx,2),zeros(numel(idx),1),bxPk(idx,2),byPk(idx,2),colors(4,:));
  
else
  % plot prominence
  plotLines(hAxes,'Prominence',x(iPk), y(iPk), x(iPk), bPk, colors(2,:));  
  
  % plot width
  plotLines(hAxes,'HalfProminenceWidth',wxPk(:,1), hm, wxPk(:,2), hm, colors(3,:));
  
  % plot peak borders
  idx = find(bPk(:)<byPk(:,1));
  plotLines(hAxes,'Border',bxPk(idx,1),bPk(idx),bxPk(idx,1),byPk(idx,1),colors(4,:));
  idx = find(bPk(:)<byPk(:,2));
  plotLines(hAxes,'Border',bxPk(idx,2),bPk(idx),bxPk(idx,2),byPk(idx,2),colors(4,:));
end

if coder.target('MATLAB')
  hLine = get(hAxes,'Children');
  tags = get(hLine,'tag');
  
  legendStrs = {};
  searchTags = {'Signal','Peak','Prominence','Height','HalfProminenceWidth','HalfHeightWidth','Border'};
  for i=1:numel(searchTags)
    if any(strcmp(searchTags{i},tags))
      legendStrs = [legendStrs, ...
        {getString(message(['signal:findpeaks:Legend' searchTags{i}]))}]; %#ok<AGROW>
    end
  end
  
  if numel(hLine)==1
    legend(getString(message('signal:findpeaks:LegendSignalNoPeaks')), ...
      'Location','best');
  else
    legend(legendStrs,'Location','best');
  end  
end

end
%--------------------------------------------------------------------------
function plotLines(hAxes,tag,x1,y1,x2,y2,c)
% concatenate multiple lines into a single line and fencepost with NaN
n = numel(x1);
line(reshape([x1(:).'; x2(:).'; NaN(1,n)], 3*n, 1), ...
     reshape([y1(:).'; y2(:).'; NaN(1,n)], 3*n, 1), ...
     'Color',c,'Parent',hAxes,'tag',tag);

end
%--------------------------------------------------------------------------
function scalePlot(hAxes)

% In the event that the plot has integer valued y limits, 'axis auto' may
% clip the YLimits directly to the data with no margin.  We search every
% line for its max and minimum value and create a temporary annotation that
% is 10% larger than the min and max values.  We then feed this to "axis
% auto", save the y limits, set axis to "tight" then restore the y limits.
% This obviates the need to check each line for its max and minimum x
% values as well.

minVal = Inf;
maxVal = -Inf;

if coder.target('MATLAB')
  hLines = findall(hAxes,'Type','line');
  for i=1:numel(hLines)
    data = get(hLines(i),'YData');
    data = data(isfinite(data));
    if ~isempty(data)
      minVal = min(minVal, min(data(:)));
      maxVal = max(maxVal, max(data(:)));
    end
  end
  
  axis auto
  xlimits = xlim;
  
  % grow upper and lower y extent by 5% (a total of 10%)
  p = .05;    
  y1 = (1+p)*maxVal - p*minVal;
  y2 = (1+p)*minVal - p*maxVal;
  
  % artificially expand the data range by the specified amount
  hTempLine = line(xlimits([1 1]),[y1 y2],'Parent',hAxes);  
  
  % save the limits
  ylimits = ylim;
  delete(hTempLine);  
else
  axis auto
  ylimits = ylim;
end

% preserve expanded y limits but tighten x axis.
axis tight
ylim(ylimits);  

% [EOF]

end


function chkinputdatatype_In(varargin)
%CHKINPUTDATATYPE Check that all inputs are double

%   Copyright 2009-2013 The MathWorks, Inc.
for n = 1:nargin
    if ~isa(varargin{n},'double')
        error(message('signal:chkinputdatatype:NotSupported'));
    end
end
end

function states = bz_LoadStates_StateEditorWrapper_In(basePath,timevector)
% TheStateEditor-appropriate loading function for loading buzcode SleepState.states.mat files.
% Abstracted here so it can be used both during initial load and during LoadStates calls
% Brendon Watson 4/2018
baseName=bz_BasenameFromBasepath(basePath);

SleepState = bz_LoadStates(basePath,'SleepState');
if isfield(SleepState,'idx')
    %Interpolate to the StateEditor timestamps
    states = interp1(SleepState.idx.timestamps,SleepState.idx.states,timevector,'nearest');
elseif isfield(SleepState,'ints')
    %If no idx saved... get from  .ints and save with idx
    SleepState.idx = bz_INTtoIDX(SleepState.ints,'statenames',{'WAKE','','NREM','','REM'});
    save(fullfile(basePath,[baseName '.SleepState.states.mat']),'SleepState') %Save with new idx
    %Interpolate to the StateEditor timestamps
    states = interp1(SleepState.idx.timestamps,SleepState.idx.states,timevector,'nearest');
elseif isempty(SleepState)
    states = zeros(1,length(timevector));
else
   error('Your SleepState is broken.')
end

states(isnan(states)) = 0; %Bug fix DL

%% END for TheStateEditor
end
