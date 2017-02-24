function [detected_swr] = detect_swr(Filebase, Channels, Epochs, varargin)
% DETECT_SWR    Detect Sharp Wave Ripple in epochs of continuous data. This
%               code is designed for linear probes, with at least one shank
%               that simultaneously records across the pyramidal cell layer
%               i.e. the superficial contacts record the ripple and at
%               least the deepest contact records the Sharp Wave (see the
%               .png files within the detect_swr for examples).
%                   |o  |
%                   |  o| <-- superficial channels record ripple
%                   |o  |
%                   |  o|
%                   |o  |
%                   |  o|
%                   |o  |
%                   |  o|  <-- deep channels record Sharp Wave
%                   \   /
%                    \ /
%                     .
%               Before running this reocord inspect your .eeg or .lfp file
%               to find that shank which best shows the sharp wave/ripple
%               profile across the shank e.g. Channels = 25:32.
%
% USE CASES
% 1) uses default parameters
% detected_swr = detect_swr( Filebase, Channels, Epochs )
%
% 2) supply parameters as a structure with fields an values
% detected_swr = detect_swr( ..., options )
%
% 3) supply parameters as Property/Value pairs
% detected_swr = detect_swr( ..., 'swBP', [2 50], 'WinSize', 300, etc. )
%   Filtering routines adapted from 'detect_hfos' by Eran Stark
%
% NOTE: -see MEMORY CHECK for how this routine handles potential memory
%       overflows associated with reading in large field potential files
%       a.k.a how I try to avoid running over into swap files
%       
%       -If you would like to see how the detector performs against a hand
%       labeled .evt file you created: name the .evt file Filebase.rip.evt 
%       with events 'start' and 'stop'.
%
%       -To improve detection, statistical thresholds are set using a local
%       window around each putative detection. This means global changes in
%       SWR and ripple magnitude aren't documented. Though, in exchange, this
%       routine is robust to non-stationarities in the data.
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% Filebase:     can be: -full or relative path with or without .xml
%                       extension
%               e.g. /mnt/Data/buddy140_060813_reo/buddy140_060813_reo.xml
%                    buddy140_060813_reo.xml
%                    buddy140_060813_reo
%
% Channels:     an array of channel ids following LoadBinary format(base 1)
%                   The array of Channels is read from superficial to deep
%                   left to right i.e. the numbering can be random,
%                   reflecting your anatomical group numbering e.g. [2 4 1]
%
% Epochs:       a list of time stamps [start stop] in seconds demarcating 
%               epochs in which to detect SWR. If this argument is empty, 
%               the entire .eeg/.lfp file will be used.               
%
%%%%%%%%%%%%%%%
%%% OUTPUTS %%%
%%%%%%%%%%%%%%%
% detected_swr: structure with fields
%
%            Ts: SWR events as [peak start stop] in seconds
%         SwMax: Sharp-wave magnitude [zscore "local empirical percentile"]
%        RipMax: Ripple magnitude [zscore "local empirical percentile"]
%        params: a structure with fields corresponding to all the
%               parameters used in the analysis
%
% OPTIONAL
% event file:   an event file formatted as filename.R01.evt that is
%               compatible with neuroscope for inspecting detections 
%               [peak start stop]
% 
% The analysis also writes out log file (.txt) containing the code used.
%
% DEPENDENCIES: all the major routines are contained within the associated
%               private folder. All should need to do is add the detect_swr
%               folder to your path. One of the training processes runs 
%               linear discriminant analysis on the joint and marginal 
%               spaces of the features. There is a statistics toolbox 
%               dependency, R2013a, which I believe requires different 
%               commands in later versions. I can update this if this 
%               routine becomes widely used.
%
%%%%%%%%%%%%%%%%%%%%%%%
%%% FREE PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%
% swBP:         [low high], passband for filtering sharp wave activity
% ripBP:        [low high] passband for filtering around ripple activity
% per_thresswD: a threshold placed upon the sharp wave difference magnitude
%               of the candidate SWR cluster determined via k-means.
% per_thresRip: a threshold placed upon ripple power based upon the non-SWR
%               cluster determined via k-means.
% WinSize:      window size in milliseconds for non-overlapping detections
% Ns_chk:       sets a window [-Ns_chk, Ns_chk] in seconds around a
%               candidate SWR detection for estimating local statistics.
% thresSDswD:   [low high], high is a threshold in standard deviations upon
%               the detected maximum sharp wave difference magnitude, based
%               upon the local distribution. low is the cutoff for the 
%               feature around the detected peak for determining the 
%               duration of the event, also based upon the local
%               distribution.
% thresSDrip:   [low high], high is a threshold in standard deviations upon
%               the detected maximum ripple magnitude, based
%               upon the local distribution. low is the cutoff for the 
%               feature around the detected peak for determining the 
%               duration of the event, also based upon the local
%               distribution.
% minIsi:       a threshold setting the minimum time between detections in
%               seconds.
% minDurSW:     a threshold setting the minimum duration of a localized
%               sharp-wave event.
% maxDurSW:     a threshold setting the maximum duration of a localized
%               sharp-wave event.
% minDurRP:     a threshold setting the minimum duration of a localized
%               ripple event associated with a sharp-wave.
% EVENTFILE:    boolean, a flag that triggers whether or not to write out
%               an event file corresponding to the detections (for
%               inspection in neuroscope).
% FIGS:         boolean, a flag that allows the user to output a figure
%               detailing the final clustering of the features and SWR
% TRAINING:     boolean, a flag for processing labeled training data
%               supplied by the user.
% DEBUG:        boolean, a flag that triggers the program to output an
%               additional set of figures and information.
% MAXGIGS:      the maximum number of GIGABYTES this 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET DEFAULT FREE PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defaults.swBP           = [ 2 50 ];
Defaults.ripBP          = [ 80 250 ];
Defaults.per_thresswD   = 10;         % percentile
Defaults.per_thresRip   = 50;         % percentile
Defaults.WinSize        = 200;        % ms
Defaults.Ns_chk         = 5;          % seconds
Defaults.thresSDswD     = [0.5 2.5];  % standard deviations
Defaults.thresSDrip     = [0.5 2.5];  % standard deviations
Defaults.minIsi         = 0.10;       % seconds
Defaults.minDurSW       = 0.02;       % seconds
Defaults.maxDurSW       = 0.5;        % seconds
Defaults.minDurRP       = 0.025;      % seconds
Defaults.EVENTFILE      = true;
Defaults.FIGS           = false;
Defaults.TRAINING       = false;
Defaults.DEBUG          = false;
Defaults.MAXGIGS        = 12;        % GigaBytes
%%%%%%%%%%%%%%%%%%%%
%%% INPUT CHECKS %%%
%%%%%%%%%%%%%%%%%%%%
mfname = mfilename;
% 1) fileabase
% user might have input an absolute or relative path
[pathname, filename, extname] = fileparts(Filebase);
if isempty(pathname)
    pathname = pwd;
end

% check whether user provided a file extension
if isempty(extname)
    filechk = dir([Filebase '.xml']);
    if isempty(filechk)
        error(['%s: incompatible format for 1st argument (Filebase)\n', ...
               '\tIf user supplied a relative path, the .xml and all\n', ...
               '\trelevant files must be in the current directory.\n'],mfname);
    end
else
    filechk = dir(Filebase);
    if isempty(filechk)
        error(['%s: incompatible format for 1st argument (Filebase)\n', ...
               '\tIf user supplied a relative path, the .xml and all\n', ...
               '\trelevant files must be in the current directory.\n'],mfname);
    end
end
Filebase = [pathname filesep filename];

% Read in xml of file parameters
fprintf(1, 'Loading XML file...\n');
par        = LoadXml( [Filebase '.xml'] );
% pull parameters from .xml file
SR         = par.lfpSampleRate; % lfp sampling rate
nChan      = par.nChannels;     % number of channels in the recording

% 2) Channels
if ~all(Channels > 0 & Channels < nChan)
    error(['%s: incompatible input for 2nd argument (Channels)\n', ...
            '\tAccording to the .xml file, the user did not supply\n', ...
            '\ta valid list of channels.\n'],mfname);
end
Nchan       = length(Channels);

% Access lfp looking for an .eeg file first and a .lfp file second
if ~isempty(dir([Filebase '.eeg']))
    lfp_file = [Filebase '.eeg'];
    lfp_info = dir(lfp_file);
elseif ~isempty(dir([Filebase '.lfp']))
    lfp_file = [Filebase '.lfp'];
    lfp_info = dir(lfp_file);
else
    error(['%s: Field potential file could not be found\n', ...
          '\tThere is neither a .eeg nor a .lfp file matching the\n', ...
          '\tbase name of the directory supplied by the user.\n'],mfname);
end

% 3) Epochs
if nargin < 3 || isempty(Epochs)
    % use the whole recording if no restrictions supplied by user
    fprintf( 1, '%s: Using Default Epochs:\n', mfname );
    fprintf( 1, '\tNo Epochs supplied by user. Using the whole recording.\n');
    Epochs = [1/SR lfp_info.bytes/(2*nChan*SR)];
elseif ~all(Epochs(:,1) < Epochs(:,2)) || ...
        ~all(Epochs(:) > 0 & Epochs(:) < lfp_info.bytes/(2*nChan*SR)) || ...
        size(Epochs,2) > 2
    error(['%s: Epochs are not formatted properly:\n', ...
           '\tEpochs must be formatted as a list of [start stop] times in seconds.\n'], mfname);
end
Nepochs = size(Epochs,1);

%%%%%%%%%%%%%%%%%%%%%%%
%%% PARSE ARGUMENTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin)
    % There is no argument validation herein.
    Defaults = parseArgs(Defaults,varargin);
end
% return parameters stored in Defaults
v2struct(Defaults);

%%%%%%%%%%%%%%%%%%%%
%%% MEMORY CHECK %%%
%%%%%%%%%%%%%%%%%%%%
% This routine requires 4*NumberOfChannelsRequested + 2 GigaBytes.
% This estimate includes memory used to copy the inputs into functions (my
% kingdom for a pointer!).
if isfield(par,'nBits')
    ByteConversion = 64/par.nBits;  % (double in Bytes)/(nBits -> Bytes)
else
    ByteConversion = 4;
end
BytesPerChan   = (lfp_info.bytes/nChan)*ByteConversion;
GigsRequested  = (4*BytesPerChan*Nchan + 2*BytesPerChan)/(1024^3);
if GigsRequested > MAXGIGS
    error(['%s: Too much memory requested\n', ...
           '\tThe user requested %4.1f GigaBytes and the maximum\n', ...
           '\tcurrently allowed is set by the variable MAXGIGS at %4.1f GigaBytes.\n', ...
           '\tIf your system can handle more memory, increase the value of MAXGIGS.\n'],mfname,GigsRequested,MAXGIGS);
end

%%%%%%%%%%%%%%%%%%%%

if TRAINING
    % Access hand labeled ripple events and convert to seconds
    fprintf( 1, 'Loading Labeled Event File for Training...\n');
    
    training_rip  = [Filebase '.rip.evt'];
    if isempty(dir(training_rip))
        error(['%s: Training Event file could not be found\n', ...
                '\tThere is not an event file "Filebase.rip.evt" in \n', ...
                '\tthe directory supplied by the user.\n'],mfname);
    end
    ripfile       = fopen(training_rip,'r');

    % count the number of lines in the file and rewind
    count = 0;
     while ischar(fgetl(ripfile))
         count = count + 1;
     end
    frewind(ripfile);

    % Access ripple events and convert to seconds [start, stop]
    if mod(count,2)
        fclose(ripfile);
        error('Check your event file: you should have start/stop pairs, but you have an odd number of lines');
    end

    % Parse event file to pull out hand labeled start/stop event pairs
    rippleTs = zeros(count/2,2);
    count     = 1;
    while 1
        tline = fgetl(ripfile);
        if ~ischar(tline), break, end
        
        if strfind(tline,'start')

            ws                         = regexp(tline,'\s');
            rippleTs(round(count/2),1) = str2double(tline(1:ws-1))/1000;

        elseif strfind(tline,'stop')

            ws                         = regexp(tline,'\s');
            rippleTs(round(count/2),2) = str2double(tline(1:ws-1))/1000;
        end
        count = count + 1;
    end
    fclose(ripfile);

    % For each ripple, I will select a segment of lfp data of the same length,
    % that does not include a ripple. I will use these paired samples for
    % comparison.
    nonrippleTs    = zeros(size(rippleTs));
    counter        = 0;

    % for each trial, find the labeled ripples contained
    fprintf( 1, 'Parsing non-ripple paired samples...\n');
    for ii = 1:Nepochs

        ripples = find(rippleTs(:,1) > Epochs(ii,1) & rippleTs(:,1) < Epochs(ii,2));
        if isempty(ripples)
            continue
        end

        % for each ripple within this trial, select a corresponding nonripple
        % segment, checking to make sure it doesn't overlap with any ripples
        % within that trial.
        Trip           = rippleTs(ripples,:);
        Nrip           = length(ripples);
        nonripTmp      = zeros(Nrip,2);
        nonripTmp(:,1) = (Epochs(ii,2)-Epochs(ii,1))*rand(Nrip,1) + Epochs(ii,1);
        nonripTmp(:,2) = nonripTmp(:,1) + diff(Trip,[],2);

        % Check nonripple segments for overlap with ripple segments and each
        % other
        flags = true(Nrip,1);
        inds  = 1:Nrip;
        while any(flags)
            for jj = inds
                if flags(jj)

                    % Concatenate all nonripples and all ripple times, 
                    % excluding the one under inspection.
                    if Nrip > 1
                        temp               = nonripTmp;
                        temp(jj,:)         = [];    % exclude current nonripple
                        Temp               = zeros(2*Nrip - 1,2);
                        Temp(1:Nrip - 1,:) = temp;
                        Temp(Nrip:end,:)   = Trip;
                    else
                        Temp = Trip;
                    end

                    % Check for overlap between current "nonripple" and ripples
                    % and other "nonripples"
                    if isempty(find(nonripTmp(jj,1) > Temp(:,1) & nonripTmp(jj,1) < Temp(:,2),1)) && ...
                       isempty(find(nonripTmp(jj,2) > Temp(:,1) & nonripTmp(jj,2) < Temp(:,2),1))

                        % log nonripple
                        nonrippleTs(counter+inds(jj),:) = nonripTmp(jj,:);
                        flags(jj)                       = 0;

                    else

                        % Resample
                        nonripTmp(jj,1) = (Epochs(ii,2)-Epochs(ii,1))*rand + Epochs(ii,1);
                        nonripTmp(jj,2) = nonripTmp(jj,1) + diff(Trip(jj,:));

                        % set all flags to true to double check that this
                        % resample didn't result in overlap
                        flags = true(Nrip,1);
                    end
                end
            end
        end
        % increment counter
        counter = counter + Nrip;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ESTIMATE LFP FEATURES: START %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As a consession to memory usage, I interleave loading of field potential
% variables, and extracting segments in the case of training data.

% load lfp file
fprintf( 1, 'Loading LFP file...\n');
lfp        = LoadBinary(lfp_file,'frequency',SR,'nChannels',nChan,...
          'channels',Channels);

% filter for sharp wave using Eran's routines and compute SW diff feature
hSw1       = makegausslpfir( swBP(1), SR, 6 );
hSw2       = makegausslpfir( swBP(2), SR, 6 );

lfpLow     = firfilt( lfp, hSw2 );      % lowpass filter
eegLo      = firfilt( lfpLow, hSw1 );   % highpass filter
lfpLow     = lfpLow - eegLo;            % difference of Gaussians
clear eegLo

swDiff     = lfpLow(:,1) - lfpLow(:,end);

if TRAINING
    % Extract ripple and paired non-ripple segments from lfp
    fprintf( 1,'Training: Extracting Paired Segments of lfp, lfpLow, and swDiff...\n');
    
    % convert timestamps to samples
    rippleSamp    = floor(rippleTs*SR);
    nonrippleSamp = floor(nonrippleTs*SR);

    % cut out ripple and paired non-ripple segments
    pairs0 = cell(counter,2);
    pairs1 = cell(counter,2);
    pairs2 = cell(counter,2);
    
    for ii = 1:counter
        pairs0{ii,1} = lfp(rippleSamp(ii,1):rippleSamp(ii,2),:);
        pairs0{ii,2} = lfp(nonrippleSamp(ii,1):nonrippleSamp(ii,2),:);

        pairs1{ii,1} = lfpLow(rippleSamp(ii,1):rippleSamp(ii,2),:);
        pairs1{ii,2} = lfpLow(nonrippleSamp(ii,1):nonrippleSamp(ii,2),:);

        pairs2{ii,1} = swDiff(rippleSamp(ii,1):rippleSamp(ii,2),1);
        pairs2{ii,2} = swDiff(nonrippleSamp(ii,1):nonrippleSamp(ii,2),1);
        
    end
end
clear lfpLow

% I made 2 changes to Eran's defaults for ripple processing:
% 1) I subtract the across channel mean from each channel prior to
% filtering
% 2) After subtracting, filtering, envelope, powerwin, and take the maximum
% across all channels at each time point as the ripple power.

% filter and process ripples, using Eran's defaults (diff of Gaussians)
hRip1      = makegausslpfir( ripBP( 1 ), SR, 6 );
hRip2      = makegausslpfir( ripBP( 2 ), SR, 6 );

% subtract mean before filtering
lfpM       = lfp - repmat(mean(lfp,2),[1 Nchan]);
clear lfp

rip        = firfilt( lfpM, hRip2 );    % highpass filter
clear lfpM
eegLo      = firfilt( rip, hRip1 );     % lowpass filter
rip        = rip - eegLo;               % difference of gaussians
clear eegLo

ripWindow  = pi / mean( ripBP );
powerWin   = makegausslpfir( 1 / ripWindow, SR, 6 );

rip        = abs(rip);
ripPower0  = firfilt( rip, powerWin );
ripPower0  = max(ripPower0,[],2);
clear rip

if TRAINING
    % Extract ripple and paired non-ripple segments from lfp
    fprintf( 1,'Training: Extracting Paired Segments of ripple Power...\n');
    
    % cut out ripple and paired non-ripple segments
    pairs3 = cell(counter,2);

    for ii = 1:counter
        
        pairs3{ii,1} = ripPower0(rippleSamp(ii,1):rippleSamp(ii,2),1);
        pairs3{ii,2} = ripPower0(nonrippleSamp(ii,1):nonrippleSamp(ii,2),1);

    end
    
    if DEBUG
        % Inspect paired segements
        figure('Position',get(0,'ScreenSize'))
        % Examine pairs, including low pass filtered data
        for ii = 1:counter
            set(gcf,'Name',sprintf('SWR Pair %d',ii))
            for jj = 1:Nchan
                
                h1 = subplot(Nchan + 2, 2, 2*jj - 1); ...
                    plot(pairs0{ii,1}(:,jj)), hold on
                    plot(pairs1{ii,1}(:,jj),'r'), hold off
                    ylims1 = get(h1,'ylim');
                    
                h2 = subplot(Nchan + 2, 2, 2*jj); ...
                    plot(pairs0{ii,2}(:,jj)), hold on
                    plot(pairs1{ii,2}(:,jj),'r'), hold off
                    ylims2 = get(h2,'ylim');
                
                if jj == 1
                    set(get(h1,'title'),'string','Labeled SWR')
                    set(get(h2,'title'),'string','Paired Non-SWR sample')
                end
                set(h1,'ylim',[min([ylims1(1) ylims2(1)]) max([ylims1(2) ylims2(2)])])
                set(h2,'ylim',[min([ylims1(1) ylims2(1)]) max([ylims1(2) ylims2(2)])])
            end
            
            h3 = subplot(Nchan + 2, 2, 2*Nchan + 1);
                plot(pairs2{ii,1},'b'), hold off
                ylims3 = get(h3,'ylim');
                ylabel('SW Diff')
            h4 = subplot(Nchan + 2, 2, 2*Nchan + 2);
                plot(pairs2{ii,2},'r'), hold off
                ylims4 = get(h4,'ylim');
            h5 = subplot(Nchan + 2, 2, 2*(Nchan + 1) + 1);
                plot(pairs3{ii,1},'b'), hold off
                ylims5 = get(h5,'ylim');
                ylabel('Ripple Power')
            h6 = subplot(Nchan + 2, 2, 2*(Nchan + 1) + 2);
                plot(pairs3{ii,2},'r'), hold off
                ylims6 = get(h6,'ylim');
            
            set(h3,'ylim',[min([ylims3(1) ylims4(1)]) max([ylims3(2) ylims4(2)])])
            set(h4,'ylim',[min([ylims3(1) ylims4(1)]) max([ylims3(2) ylims4(2)])])
            set(h5,'ylim',[min([ylims5(1) ylims6(1)]) max([ylims5(2) ylims6(2)])])
            set(h6,'ylim',[min([ylims5(1) ylims6(1)]) max([ylims5(2) ylims6(2)])])
            
            drawnow
            pause(.1)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ESTIMATE LFP FEATURES: END %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TRAINING
    
    fprintf(1,'Estimating Performance on Training Data...\n');
    
    % Assess how informative these features are for detection
    % as a first pass, I will take the max of each feature with each paired
    % sample window
    swDiffMax   = zeros(counter,2);
    ripPowerMax = zeros(counter,2);

    for ii = 1:counter
        swDiffMax(ii,1)   = max(pairs2{ii,1});
        swDiffMax(ii,2)   = max(pairs2{ii,2});

        ripPowerMax(ii,1) = max(pairs3{ii,1});
        ripPowerMax(ii,2) = max(pairs3{ii,2});
    end

    % Calculate marginal histograms
    [n1a, ctr1a] = hist(swDiffMax(:,1),40);
    [n1b, ctr1b] = hist(swDiffMax(:,2),40);
    [n2a, ctr2a] = hist(ripPowerMax(:,1),40);
    [n2b, ctr2b] = hist(ripPowerMax(:,2),40);

    % run a basic k-means clustering and evaluate performance
    [idx,Ctrs]  = kmeans([swDiffMax(:) ripPowerMax(:)],2);
    
    labels      = [ones(counter,1);2*ones(counter,1)];
    scores      = logical(idx == labels);

    % let's run a linear discriminant analysis on the joint and marginal spaces
    % (requires the statistics toolbox, R2013a, I believe the commands change 
    % in later versions)
    FeaturesAll = [swDiffMax(:) ripPowerMax(:)];
    cls1  = ClassificationDiscriminant.fit(FeaturesAll,labels); % joint
    cls2  = ClassificationDiscriminant.fit(FeaturesAll(:,1),labels); % swDiff
    cls3  = ClassificationDiscriminant.fit(FeaturesAll(:,2),labels); % rip Power
    % access performance on training data
    perf0 = 100 - ((resubLoss(cls1)*cls1.NObservations)/cls1.NObservations)*100;
    perf1 = 100 - ((resubLoss(cls2)*cls2.NObservations)/cls2.NObservations)*100;
    perf2 = 100 - ((resubLoss(cls3)*cls3.NObservations)/cls3.NObservations)*100;

    % For the unsupervised case, I'll set thresholds for ripple power
    % and sharp wave difference magnitude according to thresholds upon the
    % empirical distributions. This requires I map the probabilities of one
    % cluster to that of the other.
    x       = ripPowerMax(:,1);
    [Fi,xi] = ecdf(x);
    n       = length(x);
    xj      = xi(2:end);
    Fj      = (Fi(1:end-1)+Fi(2:end))/2;
    xj      = [xj(1)-Fj(1)*(xj(2)-xj(1))/((Fj(2)-Fj(1)));
               xj;
               xj(n)+(1-Fj(n))*((xj(n)-xj(n-1))/(Fj(n)-Fj(n-1)))];
    Fj   = [0; Fj; 1];
    
    if DEBUG
        figure('Name','Empirical Cumulative Probability Distribution')
        stairs(xi,Fi,'r');
        hold on
        plot(xj,Fj,'b-');
        hold off
    end

    % Map the Nth percentile in the non-SWR cluster to the Pth in the SWR
    % cluster.
    per_thres = 50;
    per_corr  = 100*interp1(xj,Fj,prctile(ripPowerMax(:,2),per_thres));

    % Main Figure in the case of hand labeled training data
    % close all
    figure('Color','w'),

    string = sprintf(['Sharp Wave Difference vs. Ripple Power\n',...
              'Mean subtraction before filter\n', ...
              'Max at each time point is taken for ripple power']);
    
    % Joint space of SWR difference and ripple power
    h1 = subplot(2,2,2);
        plot(swDiffMax(:,1), ripPowerMax(:,1),'.b'), hold on
        plot(swDiffMax(:,2), ripPowerMax(:,2),'.r'), hold on
        axis([min(swDiffMax(:)) max(swDiffMax(:)) 0 max(ripPowerMax(:))])
        % LDA coefficients
        K1 = cls1.Coeffs(1,2).Const;
        L1 = cls1.Coeffs(1,2).Linear;
        % LDA function
        f1 = @(x1,x2) K1 + L1(1)*x1 + L1(2)*x2;
        % Percentile threshold
        line(get(h1,'xlim'), repmat(prctile(ripPowerMax(:,1),per_corr), [1 2]), ...
            'linestyle','--','color','k')
        l0 = ezplot(h1,f1,[get(h1,'xlim') get(h1,'ylim')]);
        set(l0, 'Color','g','linewidth',2)
        title(h1,string,'fontsize',12)
        legend('SWR','~SWR','Location','NorthEastOutside')
        xlabel('Sharp Wave Diff','fontsize',12)
        ylabel('Ripple Power','fontsize',12)
        set(h1,'fontsize',12)
    % Marginal Space of Sharp Wave Difference
    h2 = subplot(2,2,4);
        bar(ctr1a,-n1a,1,'b'), hold on, bar(ctr1b,-n1b,1,'r'), ...
        axis([min(swDiffMax(:)) max(swDiffMax(:)) -max([n1a n1b])*1.1 0]),
        % LDA coefficients
        K2 = cls2.Coeffs(1,2).Const;
        L2 = cls2.Coeffs(1,2).Linear;
        % LDA function
        f2 = @(x1,x2) K2 + L2(1)*x1;
        l1 = ezplot(h2,f2,[get(h2,'xlim') get(h2,'ylim')]);
        set(l1, 'Color','g','linewidth',2)
        title('')
        xlabel('')
        ylabel('')
        axis('off'),
    % Marginal Space of Ripple Power
    h3 = subplot(2,2,1);
        barh(ctr2a,-n2a,1,'b'), hold on, barh(ctr2b,-n2b,1,'r'), ...
        axis([-max([n2a n2b])*1.1 0 min(ripPowerMax(:)) max(ripPowerMax(:))])
        % LDA coefficients
        K3 = cls3.Coeffs(1,2).Const;
        L3 = cls3.Coeffs(1,2).Linear;
        % LDA function
        f3 = @(x2,x1) K3 + L3(1)*x1;    % note the flip in x and y
        l2 = ezplot(h3,f3,[get(h3,'xlim') get(h3,'ylim')]);
        set(l2, 'Color','g','linewidth',2)
        title('')
        xlabel('')
        ylabel('')
        axis('off'),
    set(h1,'Position',[0.25 0.35 0.55 0.55]);
    set(h2,'Position',[.25 .1 .55 .15]);
    set(h3,'Position',[.05 .35 .15 .55]);

    string = sprintf(['K-means Clustering:\n', ...
                      '    Center1: X  = %7.2f, Y = %7.2f\n', ...
                      '    Center2: X  = %7.2f, Y = %7.2f\n', ...
                      '    Performance = %4.1f percent\n\n', ...
                      'LDA, joint space:\n',...
                      '    Performance = %4.1f percent\n\n', ...
                      'LDA, 1D SW Diff:\n',...
                      '    Performance = %4.1f percent\n\n', ...
                      'LDA, 1D Ripple Power:\n',...
                      '    Performance = %4.1f percent\n\n', ...
                      '%dth percentile of ~SWR corresponds to\n%3.1fth percentile of SWR'], ...
                      Ctrs(1,1),Ctrs(1,2),Ctrs(2,1),Ctrs(2,2), ...
                      (sum(scores)/numel(scores))*100,perf0,perf1,perf2,per_thres,per_corr);

    annotation('textbox',[0.81 0.70 .1 .1],'String',string, ...
               'edgecolor','none','fontsize',12)

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DETECT SWR EVENTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% I take the max of the sharp waves in non-overlapping windows. Once a
% maximum in the sharp wave is detected, the maximum ripple power around
% that detected is logged.

% Specifically, maximum detection is with reference to the Sharp Wave. If a
% 3 point maximum is detected for the sharp wave feature, then a maximum 
% for the ripple power is also logged for that detection within a window 
% around sharp wave detection.

% This first pass is permissive: it is just detecting maxima without 
% regards for magnitude

WinSize     = floor((WinSize*SR)/1000); %#ok<NODEF> % window in samples
HalfWinSize = floor(WinSize/2); % for logging ripple power around SW

% for each trial, block features into windows
% first pass, just count blocks for pre-allocating memory
Nblocks = 0;
for ii = 1:Nepochs
    epoch_inds = floor(Epochs(ii,:)*SR);
    blocks     = epoch_inds(1):WinSize:epoch_inds(2);
    Nblocks    = Nblocks + numel(blocks)-1;
end

% Second pass populates values
swDiffAll   = zeros(Nblocks,1);
ripPowerAll = zeros(Nblocks,1);

% In samples at SR
featureTs   = zeros(Nblocks,1);

MaxCounter  = 0;
for ii = 1:Nepochs
    epoch_inds = floor(Epochs(ii,:)*SR);
    blocks     = epoch_inds(1):WinSize:epoch_inds(2);
    
    for jj = 1:length(blocks)-1
        %Find the maximum sharp-wave power in the block
        [tempswDiff,ind0] = max(swDiff(blocks(jj):blocks(jj+1)-1,1));
        
        % Transform from sample within block to sample within record
        maxSample = blocks(jj) + ind0 - 1;
        
        % If the maximum is detected on the edges of the window, run a 3
        % point max check over the boundaries to resolve edge effect
        if (ind0 == 1 || ind0 == WinSize)
            
            % Check boundaries and then check for 3 point max
            if (maxSample == 1) || (maxSample == size(swDiff,1))
                % Don't use the boundaries of the recording
                continue
            else
                % Check for 3 point maximum
                [~, maxchk] = max([swDiff(maxSample - 1,1) ...
                                   swDiff(maxSample,1) ...
                                   swDiff(maxSample + 1,1)]);
                if abs(2 - maxchk) > 0
                    % This is not a maximum
                    continue
                end
            end
        end
        % Log features
        MaxCounter              = MaxCounter + 1;
        featureTs(MaxCounter)   = maxSample;
        swDiffAll(MaxCounter)   = tempswDiff;
        ripPowerAll(MaxCounter) = max(ripPower0(maxSample-HalfWinSize:maxSample+HalfWinSize,1));
    end
end
featureTs   = featureTs(1:MaxCounter);
swDiffAll   = swDiffAll(1:MaxCounter);
ripPowerAll = ripPowerAll(1:MaxCounter);

% run a basic k-means clustering to partition SWR from ~SWR
[idx,Ctrs] = kmeans([swDiffAll(:) ripPowerAll(:)],2);

% we know there will be far fewer SWR than not, so let's assign idx1 to
% the cluster with the fewest members i.e. SWR
if sum(idx == 1) > sum(idx == 2)
    idx1      = logical(idx == 1);
    idx2      = logical(idx == 2);
    % reverse labels
    temp      = idx1;
    idx1      = idx2;
    idx2      = temp;
    clear temp
else
    idx1      = logical(idx == 1);
    idx2      = logical(idx == 2);
end

% For the unsupervised case, I'll set thresholds for ripple power
% and sharp wave difference magnitude according to thresholds upon the
% empirical distributions.

if TRAINING
    fprintf(1,'Estimating Precision and Recall surfaces for features...');
    % Construct Precision/Recall surfaces
    thresL_swDs = 0:5:100;
    thresL_Rips = 0:5:100;

    % Calculate precision/recall values if labeled data is supplied
    % precision: true_positive / (true_positive + false_positive)
    % recall:    true_positive / (true_positive + false_negative)
    Precision = zeros([length(thresL_swDs) length(thresL_Rips)]);
    Recall    = zeros([length(thresL_swDs) length(thresL_Rips)]);

    for jj = 1:length(thresL_swDs)
        for kk = 1:length(thresL_Rips)

            per_thresswDt = thresL_swDs(jj);
            per_thresRipt = thresL_Rips(kk);

            % Sharp Wave Difference: set a lower cutoff in terms of percentiles
            thresL_swD    = prctile(swDiffAll(idx1),per_thresswDt);

            % Map the Nth percentile in the non-SWR cluster to the Pth in the SWR
            % cluster.
            thresL_Rip    = prctile(ripPowerAll(idx2,1),per_thresRipt);

            % retain SWR candidates according to clustering and thresholds
            SWRind = idx1 & swDiffAll > thresL_swD & ripPowerAll > thresL_Rip;

            % convert featureTs to seconds
            featureTsS  = featureTs*1/SR;
            SWRs        = featureTsS(SWRind);

            TP = 0; % true positive
            FP = 0; % false positive

            hits = false(size(rippleTs,1),1);
            for ii = 1:length(SWRs)

                chk = SWRs(ii) > rippleTs(:,1) & SWRs(ii) < rippleTs(:,2);

                if sum(chk) < 1
                    FP = FP + 1;
                else
                    % exclude duplicates
                    if any(hits & chk)
                        continue
                    end
                    TP = TP + 1;
                    hits = hits | chk;
                end 

            end

            FN = sum(~hits);  % false negatives

            Precision(jj,kk) = TP / (TP + FP);
            Recall(jj,kk)    = TP / (TP + FN);

        end
    end

    figure,
        subplot(121),imagesc(Precision),
            xlabel('thresL rip'),ylabel('thresL swD'),title('Precision'),
            set(gca,'ytick',1:2:21)
            set(gca,'yticklabel',0:10:100)
            set(gca,'xtick',1:2:21)
            set(gca,'xticklabel',0:10:100)
            colorbar
            axis square

        subplot(122),imagesc(Recall),
            xlabel('thresL rip'),ylabel('thresL swD'),title('Recall'),
            set(gca,'ytick',1:2:21)
            set(gca,'yticklabel',0:10:100)
            set(gca,'xtick',1:2:21)
            set(gca,'xticklabel',0:10:100)
            colorbar
            axis square
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOCALIZING EVENTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Localize events and winnow detections using the following criteria:
fprintf('Localizing and Screening Detections ...\n');
% 1) detected sharp wave difference magnitude must be > thresSDswD(2)
% 2) detected ripple power must be > thresSDrip(2)
% 3) detections can't be within a minimum inter-event threshold
% 4) both the detected sharp-wave and ripple power increase mubt persistent
    % for a minimum duration.
    minDurSWs = floor(minDurSW*SR);    % convert to samples
    minDurRPs = floor(minDurRP*SR);
% 5) The sharp-wave can not exceed a maximum duration
    maxDurSWs = floor(maxDurSW*SR);
    
% Sharp Wave Difference: set a lower cutoff in terms of percentiles
thresL_swD   = prctile(swDiffAll(idx1),per_thresswD);

% Set the value in the the non-SWR cluster at the Nth percentile as the
% lower threshold for ripple detection
thresL_Rip   = prctile(ripPowerAll(idx2,1),per_thresRip);

% retain SWR candidates according to clustering and thresholds
SWRind = idx1 & swDiffAll > thresL_swD & ripPowerAll > thresL_Rip;
SWRid  = find(SWRind);

% Access time-stamps of SWR features
SWRs        = featureTs(SWRind);
Nswr        = length(SWRs);

% Access the subset of detections in samples to check for localization
bound       = Ns_chk*SR;    % convert to samples

% Exclude candidate SWRs that don't have a +/- bounds around them. Of
% course, if you make bounds huge, or have a short recording, this will be
% a problem
exclusions = (SWRs - bound < 1) | (SWRs + bound > size(swDiff,1));
SWRs(exclusions)   = [];
SWRind(exclusions) = [];
SWRid(exclusions)  = [];

fprintf(1,'Number of candidate SWRs excluded due to boundary constraints: %d\n',Nswr - length(SWRs));
Nswr        =  length(SWRs);
timeline    = -bound:bound;
SWR_diff    =  [0; diff(SWRs*(1/SR))];    % in seconds

% localize events
SWR_valid.Ts     = zeros(size(SWRs,1),3);   % peak, start, stop in samples
SWR_valid.SwMax  = zeros(size(SWRs,1),2);   % z-score & percentile
SWR_valid.RipMax = zeros(size(SWRs,1),2);   % z-score & percentile

valid            = 0;
SWRindF          = SWRind;

if DEBUG
    figure
end

for ii = 1:Nswr
    
    if DEBUG
        set(gcf,'Name',sprintf('Detection # %d out of %d',ii,Nswr))
    end
    
    % threshold feature distributions to localize about detection
    
    % 1) detected sharp wave difference magnitude must be > thresSD locally
    swDiff_chk   = swDiff(SWRs(ii)-bound:SWRs(ii)+bound);
        MswDiff  = median(swDiff_chk);
        SDswDiff = std(swDiff_chk);
        
    if swDiff_chk(bound+1) < MswDiff + thresSDswD(2)*SDswDiff
        if DEBUG
            disp('Sharp Wave Difference magnitude too small locally')
        end
        SWRindF(SWRid(ii)) = false;
        continue
    end
    
    % 2) detected ripple power must be > thresSD locally
    ripPow_chk    = ripPower0(SWRs(ii)-bound:SWRs(ii)+bound);
        MripPow   = median(ripPow_chk);
        SDripPow  = std(ripPow_chk);
        
    if ripPowerAll(SWRid(ii)) < MripPow + thresSDrip(2)*SDripPow
        if DEBUG
            disp('Ripple Power too low locally')
        end
        SWRindF(SWRid(ii)) = false;
        continue
    end
    
    % 3) detections can't be within a minimum inter-event interval thresh
    if ii < Nswr
        if SWR_diff(ii) < minIsi
            if DEBUG
                disp('IEI is too small')
            end
            SWRindF(SWRid(ii)) = false;
            continue
        end
    end
    
    % Estimate duration of detected Sharp-Wave event
    startSWR = find(swDiff_chk(1:bound + 1) < MswDiff + thresSDswD(1)*SDswDiff,1,'last');
    stopSWR  = find(swDiff_chk(bound + 1: end) < MswDiff + thresSDswD(1)*SDswDiff,1,'first');
    durSWR   = (bound + 1 - startSWR) + stopSWR;
    
    % Estimate duration of detected ripple power increase event
    % find index of maximum ripple power increase around Sharp-Wave
    [~,maxRpInd]  = max(ripPow_chk(bound+1-HalfWinSize:bound+1+HalfWinSize));
    startRP       = find(ripPow_chk(1:bound - HalfWinSize + maxRpInd) ...
                         < MripPow + thresSDrip(1)*SDripPow,1,'last');
    stopRP        = find(ripPow_chk(bound - HalfWinSize + maxRpInd: end) ...
                         < MripPow + thresSDrip(1)*SDripPow,1,'first');
    durRP         = (bound - HalfWinSize + maxRpInd - startRP) + stopRP;
    
    % 4) events associated with detections must exceed a minimum duration
    
    if durSWR < minDurSWs
        if DEBUG
            disp('Sharp-Wave Event Duration is too short.')
        end
        SWRindF(SWRid(ii)) = false;
        continue
    end
    
    if durRP < minDurRPs
        if DEBUG
            disp('Ripple Event Duration is too short.')
        end
        SWRindF(SWRid(ii)) = false;
        continue
    end
    
    % 5) sharp-wave must not exceed a maximum duration
    if durSWR > maxDurSWs
        if DEBUG
            disp('Sharp-Wave Event Duration is too long.')
        end
        SWRindF(SWRid(ii)) = false;
        continue
    end
    
    % We have a valid detection
    valid                     = valid + 1;
    SWR_valid.Ts(valid,:)     = [SWRs(ii) (SWRs(ii)-(bound+1-startSWR)) (SWRs(ii)+stopSWR-1)];
    
    swMax                     = swDiffAll(SWRid(ii));
    ripMax                    = ripPowerAll(SWRid(ii));
    
    SWR_valid.SwMax(valid,:)  = ...
        [(swMax - MswDiff)/SDswDiff sum(swDiff_chk < swMax)/(2*bound + 1)];
    
    SWR_valid.RipMax(valid,:) = ...
        [(ripMax - MripPow)/SDripPow sum(ripPow_chk < ripMax)/(2*bound + 1)];
    
    if DEBUG
        h1 = subplot(2,3,[1 2]);
            plot(timeline,swDiff_chk), axis tight, hold on
            plot([timeline(startSWR) timeline(bound+stopSWR)],...
                [swDiff_chk(startSWR) swDiff_chk(bound+stopSWR)],'og')
            line([0 0],get(gca,'ylim'), ...
                'linestyle','--','color','k')
            line(get(gca,'xlim'),repmat(MswDiff,[1 2]), ...
                'linestyle','--','color','k')
            line(get(gca,'xlim'),repmat(MswDiff+thresSDswD(2)*SDswDiff,[1 2]),...
                'linestyle','--','color','r')
            line(get(gca,'xlim'),repmat(MswDiff-thresSDswD(2)*SDswDiff,[1 2]), ...
                    'linestyle','--','color','r')
            xlabel('Samples'),
            ylabel('A.U.'),
            title(sprintf(['Sharp Wave Difference Magnitude\n',...
                           'Sharp Wave Duration:%5.1f'],durSWR*(1000/SR))), hold off
        subplot(2,3,3), hist(swDiff_chk,100), ...
            line(repmat(MswDiff+thresSDswD(2)*SDswDiff, [1 2]),get(gca,'ylim'), ...
                'linestyle','--','color','k')
            line(repmat(swDiff_chk(bound+1), [1 2]),get(gca,'ylim'), ...
                'linestyle','--','color','r')
            xlabel('A.U.'),
            ylabel('Counts')
            title(sprintf(['Histogram of Sharp Wave Difference Magnitude\n', ...
                '%3.1f SD threshold (black) and detected value (red)'],thresSDswD(2)))
            hold off

        h2 = subplot(2,3,[4 5]);
            plot(timeline,ripPow_chk), axis tight, hold on
            plot([timeline(startSWR) timeline(bound+stopSWR)],...
                    [ripPow_chk(startSWR) ripPow_chk(bound+stopSWR)],'og')
            line([0 0],get(gca,'ylim'), ...
                'linestyle','--','color','k')
            line(get(gca,'xlim'),repmat(MripPow,[1 2]), ...
                'linestyle','--','color','k')
            line(get(gca,'xlim'),repmat(MripPow+thresSDrip(2)*SDripPow,[1 2]), ...
                'linestyle','--','color','r')
            line(get(gca,'xlim'),repmat(MripPow-thresSDrip(2)*SDripPow,[1 2]), ...
                'linestyle','--','color','r')
            xlabel('Samples'),ylabel('A.U.'),title('Rip Pow'), hold off
        subplot(2,3,6), hist(ripPow_chk,100), ...
            line(repmat(MripPow+thresSDrip(2)*SDripPow, [1 2]),get(gca,'ylim'), ...
                'linestyle','--','color','k')
            line(repmat(ripPow_chk(bound+1), [1 2]),get(gca,'ylim'), ...
                'linestyle','--','color','r')
            xlabel('A.U.'),
            ylabel('Counts')
            title(sprintf(['Histogram of Ripple Power\n', ...
                '%3.1f SD threshold (black) and detected value (red)'],thresSDrip(2)))
            hold off

        linkaxes([h1 h2],'x')
        
%         keyboard
        pause(.1)
        
        clf
    end
    
end
SWR_valid.Ts     = SWR_valid.Ts(1:valid,:);
SWR_valid.SwMax  = SWR_valid.SwMax(1:valid,:);
SWR_valid.RipMax = SWR_valid.RipMax(1:valid,:);

if FIGS
    % Evaluate fixed point Precision/Recall performance
    figure('Color','w'),

    % marginal histograms
    [n1, ctr1] = hist(swDiffAll,1000);
    [n2, ctr2] = hist(ripPowerAll,1000);

    h1 = subplot(2,2,2);
        plot(swDiffAll(idx1), ripPowerAll(idx1),'.r'), hold on
        plot(swDiffAll(idx2), ripPowerAll(idx2),'.b'),
        plot(swDiffAll(SWRindF),ripPowerAll(SWRindF),'og','markerfacecolor','none')
        axis([min(swDiffAll(:)) max(swDiffAll(:)) min(ripPowerAll(:)) max(ripPowerAll(:))])
        line(get(h1,'xlim'), repmat(thresL_Rip,[1 2]),'linestyle','--','color','k')
        line(repmat(thresL_swD,[1 2]),get(h1,'ylim'),'linestyle','--','color','k')
        title('Sharp Wave Difference vs. Ripple Power','fontsize',12)
        xlabel('Sharp Wave Diff','fontsize',12)
        ylabel('Ripple Power','fontsize',12)
    h2 = subplot(2,2,4);
        bar(ctr1,-n1,1),
        axis([min(swDiffAll(:)) max(swDiffAll(:)) -max(n1)*1.1 0]),
        axis('off'),
    h3 = subplot(2,2,1);
        barh(ctr2,-n2,1),
        axis([-max(n2)*1.1 0 min(ripPowerAll(:)) max(ripPowerAll(:))])
        axis('off'),
    set(h1,'Position',[0.25 0.35 0.55 0.55],'fontsize',12);
    set(h2,'Position',[.25 .1 .55 .15],'fontsize',12);
    set(h3,'Position',[.05 .35 .15 .55],'fontsize',12);

    string = sprintf(['K-means Clustering:\n\n',...
                      '    Center1: X  = %7.2f, Y = %7.2f\n', ...
                      '    Center2: X  = %7.2f, Y = %7.2f\n\n', ...
                      'Number of Detected SWR:%5d\n'],...
                      Ctrs(1,1),Ctrs(1,2),Ctrs(2,1),Ctrs(2,2), valid);
    annotation('textbox',[0.81 0.70 .1 .1],'String',string, ...
               'edgecolor','none','fontsize',12)
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATE OUTPUT %%%
%%%%%%%%%%%%%%%%%%%%%%%
% 1) swr detections
detected_swr.Ts     = SWR_valid.Ts/SR;
detected_swr.SwMax  = SWR_valid.SwMax;
detected_swr.RipMax = SWR_valid.RipMax;

%Write event file: Generate Output and Write out Event file
if EVENTFILE
    rippleFiles = dir('*.R*.evt');
    if isempty(rippleFiles)
        fileN = 1;
    else
        %set file index to next available value\
        pat = '.R[0-9].';
        fileN = 0;
        for ii = 1:length(rippleFiles)
            token  = regexp(rippleFiles(ii).name,pat);
            val    = str2double(rippleFiles(ii).name(token+2:token+4));
            fileN  = max([fileN val]);
        end
        fileN = fileN + 1;
    end
    fid = fopen(sprintf('%s%s%s.R%02d.evt',pathname,filesep,filename,fileN),'w');

    % convert detections to milliseconds
    SWR_final   = SWR_valid.Ts*(1000/SR);
    fprintf(1,'Writing event file ...\n');
    for ii = 1:size(SWR_final,1)
        fprintf(fid,'%9.1f\tstart\n',SWR_final(ii,2));
        fprintf(fid,'%9.1f\tpeak\n',SWR_final(ii,1));
        fprintf(fid,'%9.1f\tstop\n',SWR_final(ii,3));
    end
    fclose(fid);
end

% 3) analysis params
params.filename     = filename;
params.Channels     = Channels;
params.SR           = SR;
params.nChan        = nChan;
params.swBP         = swBP;
params.ripBP        = ripBP;
params.per_thresswD = per_thresswD;
params.per_thresRip = per_thresRip;
params.WinSize      = WinSize;
params.Ns_chk       = Ns_chk;
params.thresSDswD   = thresSDswD;
params.thresSDrip   = thresSDrip;
params.minIsi       = minIsi;
params.minDurSW     = minDurSW;
params.maxDurSW     = maxDurSW;
params.minDurRP     = minDurRP;
params.EVENTFILE    = EVENTFILE;
params.TRAINING     = TRAINING;
params.DEBUG        = DEBUG;

detected_swr.params = params;

% save output to Filebase directory
save(sprintf('%s_detected_swr%d.mat',Filebase,fileN),'detected_swr')

% write out log file
logAnalysisFile(mfilename('fullpath'),pathname);

end %EOF