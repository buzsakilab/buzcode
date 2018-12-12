%% BUZCODE TUTORIAL: SPW-R (HPC), UP/DOWN, and State Event Detection 
 
% SUMMARY: 
% In this tutorial we outline the use of the buzcode repository to load a 
% .lfp file and detect i) hippocampal sharpwave-ripples and ii) cortical UP/DOWN 
% states. We then create a .evt file that can be read into neuroscope for 
% simultaneous viewing with .lfp data.

% ASSUMPTIONS:
% 1.   You have gone through the I/O Tutorial and have:
% i.   Downloaded git: https://github.com/
% ii.  Cloned and compiled the buzcode repo: https://github.com/buzsakilab/buzcode
% iii. Downloaded MATLAB and added buzcode to your path
% iv.  Downloaded neuroscope to view .lfp data: http://neurosuite.sourceforge.net/ 
% 2.   You are currently in the folder < /buzcode/tutorials >

% TUTORIAL OVERVIEW:
% 1.    Load Data   
% 1i.   Create/load session info
% 1iii. Convert .dat -> .lfp (downsample data from 20kHz to 1250kHz)
% 1iv.  Load .lfp file
% 2.    Ripple detection and stats
% 3.    Intro spectral processing: filters
% 4.    Plots
% 5.    Make .evt file for use with neuroscope
% 6.    UP/DOWN event detection (incomplete)
% 7.    State detection (incomplete)

% TO DO:
% Update this tutorial once have streamlined bz_Load and bz_Save
% Add a HPC-CTX recording for UP/DOWN event detection
% Add in state scoring section 

%% 1i.  Specify variables

% Confirm user is in buzcode/tutorial folder
        path = strsplit(pwd,filesep);
        if ~strcmp(path{end},'tutorials') & ~strcmp(path{end-1},'buzcode') 
            error('this script assumes you are running it from the /buzcode/tutorials folder...') 
        end

% Specify global variables       
        basePath = fullfile(cd, 'exampleDataStructs','20170505_396um_0um_merge'); % This is the path to the tutorials folder that contains all data structures associated with the given recording session.
        baseName = bz_BasenameFromBasepath(basePath); % baseName is the base file name for each of the data structures, which corresponds to the name of the folder containing the data structures (i.e. 20170505_396um_0um_merge).         
        
%% 1ii. Load sessioninfo
        sessionInfo = bz_getSessionInfo(basePath); % SessionInfo is a structure, where each field contains information about the recording session. For more information, see < https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards > . 

%% 1iv.  Load LFP data
        
% Specify hippocampal channel used for ripple detection
        try
            lfpChan = sessionInfo.ca1; % NOTE: add in fx for auto detection of good ripple channel if sessionInfo.ca1 field does not exist
        catch
            error('SessionInfo does not specify HPC channel for SPW-R detection. Manually specify by inspecting .lfp file in neuroscope -or- detect channel autmoatically using bz_FindRippleChannel')
        end
        
% Download .lfp data from online repository
        try
            disp('Loading data')
            lfp = bz_GetLFP(lfpChan,'basepath',basePath);           
        catch
            disp('Downloading data from online repository')
            cd(basePath)
            run(fullfile(basePath,'download_DATA.m')) % Download .lfp file from our online database to appropriate tutorials folder
            lfp = bz_GetLFP(lfpChan,'basepath',basePath);
        end

       % Alternative:
        % lfpAll = bz_GetLFP('all','basepath',basePath); % Load all channels
        % riplfp = lfpAll.data(:,lfp.channels==ripchan); % Example, index into structure to look at data from one channel

%% 2.   Ripple detection

% Params (defaults listed below)
        ripthresh = [2 5]; %min std, max std
        ripdur = [30 100];
        ripfreq = [130 200];
         
% Detect SPW-Rs or load structure with detection info if already detected
        if isempty(dir(fullfile(basePath, '*ripples.events.mat'))) == 1
            ripples = bz_FindRipples(basePath, lfpChan,'thresholds',ripthresh,'durations',ripdur,'show','on','saveMat',true);
        else
            ripples = bz_LoadEvents(basePath,'CA1Ripples');
        end

% Calculate ripple stats
        [maps,data,stats] = bz_RippleStats(double(lfp.data),lfp.timestamps,ripples);

% Sort by duration vs amplitude of ripple
        [~,dursort]=sort(data.duration,1,'descend');
        [~,ampsort]=sort(data.peakAmplitude,1,'descend');
        
%% 3.   Intro spectral processing: filter in ripple frequency and compare to output of bz_FindRipples

% Filter channel in ripple freq range
        ripfiltlfp = bz_Filter(lfp.data,'passband',ripfreq,'filter','fir1'); % lots of options in bz_Filter, read through it
        
%% 4.   Plots 
%% 4i.  Look: Plot some ripples
        figure()
        x=maps.ripples(dursort,:);
        for i=1:100
            subplot(10,10,i)
            plot(x(i,:))
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            axis off
        end
        
%% 4ii. Look: Detected ripples amplitude vs filtered signal (4 different views of the same thing)
        figure()
        subplot 221
        imagesc(maps.amplitude(ampsort,:))
        title('SPW-R Amplitude: sorted by amplitude')
        subplot 222
        imagesc(maps.amplitude(dursort,:))
        title('SPW-R Amplitude: sorted by duration')
        subplot 223
        imagesc(maps.ripples(ampsort,:))
        title('SPW-R Filtered Signal: sorted by amplitude')
        subplot 224
        imagesc(maps.ripples(dursort,:))
        title('SPW-R Filtered Signal: sorted by duration')
        
%% 4iii.Look: Scatterplot with marginals 
        figure()        
        scatterhist(log10(data.duration*1000),data.peakAmplitude,'kernel','on','Location','SouthWest',...
        'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
        'LineWidth',[2,2,2],'Nbins',[20 100], 'marker','.','markersize',10)
        box off
        LogScale('x',10)
        xlabel('logDuration (ms)'); ylabel('Amplitude (au)')
     
%% 5.  Make .evt file for use with neuroscope
        % This is useful for visual inspection of detected events in
        % neuroscope along with the corresponding LFP data. 
        % To load an event file, first open your data. Then, within neuroscope, go to 
        % File > Load Event File, and you will be prompted to select the .evt file.
        % Within neuroscope, a little tab will pop up beside the other two
        % in the top left area - click on that and then click on each of
        % the labels (e.g. start, peak, stop). When they are highlighted,
        % they will be displayed over your data.
       
        % Params
        eventtype = 'ripples'; %string for you that = event type you're saving
        numeventtypes = 3; % you have 3 diff types of events, start, peak, stop
        
        % Save as .evt for inspection
        ripbaseName = [baseName eventtype '.RO1.evt']; % you need the ROX bc neuroscope is buggy and uses this to parse files.
        
        % Below is for ripples specifically 
        % Populate events.time field
        lengthAll = numel(ripples.peaks)*numeventtypes;
        events.time = zeros(1,lengthAll);
        events.time(1:3:lengthAll) = ripples.timestamps(:,1);
        events.time(2:3:lengthAll) = ripples.peaks;
        events.time(3:3:lengthAll) = ripples.timestamps(:,2);
        
        % Populate events.description field
        events.description = cell(1,lengthAll);
        events.description(1:3:lengthAll) = {'start'};
        events.description(2:3:lengthAll) = {'peak'};
        events.description(3:3:lengthAll) = {'stop'};
        
        % Save .evt file for viewing in neuroscope - will save in your current directory
        SaveEvents(fullfile(basePath,ripbaseName),events) %Save and load into neuroscope along with corresponding LFP file

%% 6. UP/DOWN Event detection 
%       To do: Add cortical recording to tutorials folder for this section 

%% 7. State Scoring
%       To do: Add in this section DL
   