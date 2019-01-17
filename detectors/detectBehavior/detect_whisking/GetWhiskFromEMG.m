function [ EMGwhisk ] = GetWhiskFromEMG( basePath,varargin )
%[ EMGwhisk ] = GetWhiskFromEMG(basePath ) 
%This is a detector that extracts whisking/nonwhisking epochs from 
%implanted EMG in the whisker pad. Extracts also the EMG and EMG envelope.
%
%INPUT
%   basePath
%       Assumes presence of the following files:
%           basePath/baseName.abf    (whisking/camera pulses from clampex)
%           basePath/analogin.dat    (camera pulses to intan)
%       where basePath is a folder of the form: 
%           whateverPath/baseName/
%       If basePath not specified, tries the current path.
%
%   (parameters)
%       'gausswidth'    Gaussian smoothing width (s)        (default: 0.05)
%       'Whthreshold'   EMG Threshold for Whisking          (default: 'auto')
%       'NWhthreshold'  EMG Threshold for Non-Whisking      (default: 'auto')
%       'minwhisk'      Min. whisking duration (s)          (default: 0.1)
%       'minNWh'        Min. nonwhisking duration (s)       (default: 0.1)
%       'whiskmerge'    Min. interwhisking duration (s)     (default: 0.1)
%       'NWhmerge'      Min. internonwhisking duration (s)  (default: 0.02)
%       'showfig'       show a detection figure (default: true)
%       'saveMat'       save the output in a baseName.EMGwhisk.behavior.mat file
%       'PulseChannel'  channel with camera pulses          (default: 1)
%       'EMGChannel'    channel with EMG/Piezo signal       (default: 2)
%       'analysis'      'EMG' or 'touch'                    (default: 'EMG')
%
%   Note: 'auto' threshold detection gradient descents to a local trough in
%   the smoothed EMG envelope from an initial guess of 0.5
%
%OUTPUT
%   EMGwhisk
%   'saveMat' Creates file:
%       basePath/baseName.EMGwhisk.behavior.mat
%
%To Add: option for no camera pulses, which just align time to intan/pupil
%
%DLevenstein, WMunoz 2017
%% DEV
% basePath = '/mnt/proraidDL/Database/WMProbeData/';
% baseName = 'Layers_LFP_Test02_170323_151411';
%% Parameters
p = inputParser;
addParameter(p,'gausswidth',0.05,@isnumeric)
addParameter(p,'Whthreshold','auto',@isnumeric)
addParameter(p,'NWhthreshold','auto',@isnumeric)
addParameter(p,'minwhisk',0.1,@isnumeric)
addParameter(p,'minNWh',0.1,@isnumeric)
addParameter(p,'whiskmerge',0.1,@isnumeric)
addParameter(p,'NWhmerge',0.02,@isnumeric)
addParameter(p,'showfig',true,@islogical)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'PulseChannel',1)
addParameter(p,'EMGChannel',2)
addParameter(p,'EMGanalysis','EMG',@islogical)

parse(p,varargin{:})
SHOWFIG = p.Results.showfig;
saveMat = p.Results.saveMat;
EMGparms.gausswidth = p.Results.gausswidth;  
EMGparms.Whthreshold = p.Results.Whthreshold;
EMGparms.NWhthreshold = p.Results.NWhthreshold;
EMGparms.minwhisk = p.Results.minwhisk;
EMGparms.minNWh = p.Results.minNWh;
EMGparms.whiskmerge = p.Results.whiskmerge;
EMGparms.NWhmerge = p.Results.NWhmerge;
EMGparms.PulseChannel = p.Results.PulseChannel;
EMGparms.EMGChannel = p.Results.EMGChannel;
EMGanalysis = p.Results.EMGanalysis;

if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);
%%
abfname = fullfile(basePath,[baseName,'.abf']);
analogName = fullfile(basePath,['analogin.dat']);

switch EMGanalysis
    case 'EMG'
        savefile = fullfile(basePath,[baseName,'.EMGwhisk.behavior.mat']);
    case 'touch'
        savefile = fullfile(basePath,[baseName,'.Piezotouch.behavior.mat']);
end

figfolder = fullfile(basePath,'DetectionFigures');

if ~exist(abfname,'file')
    display('No .abf file named basePath/baseName/baseName.abf')
end
if ~exist(analogName,'file')
    display('No analogin.dat in basePath/baseName/')
end


%% Clampex File

sf_abf = 20000; %Sampling Frequency of the .abf file

[abffile,si,file_info] = abfload(abfname);

% %Prompt user for channel
% numabfchans = length(abffile(1,:));
% if ~exist('chanNums','var')
%     emgchan = listdlg('ListString',file_info.recChNames,...
%         'PromptString',['Which channels is Whisker? ']);
%     %Replace this with prompt file_info.recChNames
% end

pulse_abf = abffile(:,EMGparms.PulseChannel);
EMG = abffile(:,EMGparms.EMGChannel);

t_abf = [1:length(EMG)]'./sf_abf;

%%
% figure
% plot(t_abf,EMG,'k')
%% Get the whisking envelope
EMGrange = [400 3000];
EMG = FiltNPhase(EMG,EMGrange,sf_abf);

downsamplefactor = 16; %Downsample to same as the LFP;
EMG = downsample(EMG,downsamplefactor);
t_EMG = downsample(t_abf,downsamplefactor);
sf_down = sf_abf./downsamplefactor;

%% Z-Score the EMGZ and get EMG envelope with RMS
EMGz = NormToInt(EMG,'modZ'); %Modified Z score - robust to outliers
EMGsm = RMSEnvelope(EMGz,EMGparms.gausswidth,1/sf_down);
EMGsm = EMGsm-min(EMGsm);
EMGwhisk.EMGenvelope = EMGsm;

%% Set the thresholds by whisking troughs - ('auto')
%find by "gradient descent"(ish) from initial guess (0.5)
EMGbins = linspace(-1.5,2,100);
EMGhist = hist(log10(EMGsm),EMGbins);
EMGgrad = smooth(gradient(EMGhist),4);

%Find troughs (gradient crossing from - to +)
troughidx = find(diff(EMGgrad>0)==1);
troughs = 10.^EMGbins(troughidx);

if strcmp(EMGparms.Whthreshold,'auto')
    EMGparms.Whthreshold = 0.5;
    %Get sign of gradient at each of the thresholds and use that to pick trough
    Whsign = sign(interp1(EMGbins,EMGgrad,log10(EMGparms.Whthreshold),'nearest'));
    if Whsign==-1 
        EMGparms.Whthreshold = troughs(find(troughs>EMGparms.Whthreshold,1,'first'));
        if isempty(EMGparms.Whthreshold)
            warning('NO DIP! Try Setting manual Wh threshold, saving distribution...')
            figure
                %subplot(4,2,6)
                    hist(log10(EMGsm),100)
                    hold on
                    %plot([1 1].*log10(EMGparms.Whthreshold),get(gca,'ylim'),'g')
                    %plot([1 1].*log10(EMGparms.NWhthreshold),get(gca,'ylim'),'r')
                    axis tight
                    xlabel('EMG/Piezo Envelope (modZ)');
                    xlim([-1.5 max(log10(EMGsm))])

                    LogScale('x',10)
                    
                    switch EMGanalysis
                        case 'EMG'
                            NiceSave('FAILEDWhiskingDetection',figfolder,baseName)
                    	case 'touch'
                            NiceSave('FAILEDTouchDetection',figfolder,baseName)
                    end
            EMGwhisk = [];
            return
        end
    elseif Whsign==1
        EMGparms.Whthreshold = troughs(find(troughs<EMGparms.Whthreshold,1,'last'));
    end
end

if strcmp(EMGparms.NWhthreshold,'auto')
    EMGparms.NWhthreshold = 0.5;
    NWhsign = sign(interp1(EMGbins,EMGgrad,log10(EMGparms.NWhthreshold),'nearest'));
    if NWhsign==-1 
        EMGparms.NWhthreshold = troughs(find(troughs>EMGparms.NWhthreshold,1,'first'));
        if isempty(EMGparms.NWhthreshold)
            warning('NO DIP! Try Setting manual NWh threshold, saving distribution...')
            
            figure
                %subplot(4,2,6)
                    hist(log10(EMGsm),100)
                    hold on
                    plot([1 1].*log10(EMGparms.Whthreshold),get(gca,'ylim'),'g')
                    %plot([1 1].*log10(EMGparms.NWhthreshold),get(gca,'ylim'),'r')
                    axis tight
                    xlim([-1.5 max(log10(EMGsm))])

                    xlabel('EMG/Piezo Envelope (modZ)');
                    LogScale('x',10)
                    
                    switch EMGanalysis
                        case 'EMG'
                        NiceSave('FAILEDWhiskingDetection',figfolder,baseName)
                        case 'touch'
                        NiceSave('FAILEDTouchDetection',figfolder,baseName)
                    end
            EMGwhisk = [];
            return
        end
    elseif NWhsign==1
        EMGparms.NWhthreshold = troughs(find(troughs<EMGparms.NWhthreshold,1,'last'));
    end
end
%%
% figure
% bar(EMGbins,EMGhist)
% hold on
% plot(EMGbins,EMGgrad)
% plot([1 1].*log10(EMGparms.Whthreshold),get(gca,'ylim'),'g','LineWidth',2)
% plot([1 1].*log10(EMGparms.NWhthreshold),get(gca,'ylim'),'r','LineWidth',2)
% 
% axis tight
% xlabel('EMG Envelope (modZ)');
% LogScale('x',10)
% xlim([-1.5 max(log10(EMGsm))])

%% Identify Whisking on/offsets: EMG envelope crosses threshold
wh_thresh = EMGsm > EMGparms.Whthreshold;
wh_on = find(wh_thresh(2:end)>wh_thresh(1:end-1))+1; %whisking onsets (si)
wh_off = find(wh_thresh(2:end)< wh_thresh(1:end-1))+1;%whisking offsets (si)

% If data starts/ends in the middle of an epoch, drop first/last trigger
if wh_off(1)<wh_on(1)
    wh_off = wh_off(2:end);
end
if wh_off(end) < wh_on(end)
    wh_on = wh_on(1:end-1);
end

%% Identify NonWhisking on/offsets: EMG envelope crosses below threshold
nwh_thresh = EMGsm < EMGparms.NWhthreshold;
nwh_on = find(nwh_thresh(2:end)>nwh_thresh(1:end-1))+1; %whisking onsets (si)
nwh_off = find(nwh_thresh(2:end)< nwh_thresh(1:end-1))+1;%whisking offsets (si)

% If data starts/ends in the middle of an epoch, drop first/last trigger
if nwh_off(1)<nwh_on(1)
    nwh_off = nwh_off(2:end);
end
if nwh_off(end) < nwh_on(end)
    nwh_on = nwh_on(1:end-1);
end

%% Merge and Duration Thresholds
%Convert to seconds
wh_on = wh_on./sf_down;
wh_off = wh_off./sf_down;
nwh_on = nwh_on./sf_down;
nwh_off = nwh_off./sf_down;

%Merge brief interruptions
[ NWhints ] = MergeSeparatedInts( [nwh_on,nwh_off],EMGparms.NWhmerge );
[ Whints ] = MergeSeparatedInts( [wh_on,wh_off],EMGparms.whiskmerge );

%Drop nonwhisk epochs smaller than a minimum
[nwh_on,nwh_off] = MinEpochLength(NWhints(:,1),NWhints(:,2),EMGparms.minNWh,1);
%Drop whisking epochs smaller than a minimum
[wh_on,wh_off] = MinEpochLength(Whints(:,1),Whints(:,2),EMGparms.minwhisk,1);

%% Durations
Whdur = wh_off-wh_on;
NWhdur = nwh_off-nwh_on;

numbinst = 25;
durhist.bins = linspace(-1,2,numbinst);
durhist.Wh = hist(log10(Whdur),durhist.bins);
durhist.NWh = hist(log10(NWhdur),durhist.bins);

%% Identify Pulses from Camera (in clampex)
pulsethreshold_abf =0.5;  %Adjust this later to set based on input.
pulseonsets = find(diff(pulse_abf>pulsethreshold_abf)==1);
pulset = t_abf(pulseonsets);

interpulse = diff(pulset);

%remove the first trigger, which is just camera onset... 
%Make this more rigorous later
if interpulse(1) > sum(interpulse(2:3))
    pulset(1) = []; 
end
firstpulstime_abf = pulset(1);

%% Load the analogin for the timestamps (pulses in intan)

timepulses = bz_LoadBinary(analogName,'nChannels',1,'precision','uint16');

sf_pulse = 1./20000; %Sampling Frequency of the .abf file
t_pulse = [1:length(timepulses)]'.*sf_pulse;

pulsethreshold =1e4;  %Adjust this later to set based on input.
pulseonsets = find(diff(timepulses>pulsethreshold)==1);
pulset = t_pulse(pulseonsets);

minpulsedur = 0.003; %Remove double/noise crossings
shortpulses=diff(pulset)<(minpulsedur);
pulset(shortpulses) = [];

interpulse = diff(pulset);

if interpulse(1) > sum(interpulse(2:3))
    pulset(1) = []; 
end
sf_eff = 1./mean(interpulse);

%Check that frame duration is constant up to tolerance (no skipped frames)
tol = 0.001;

if range(interpulse)>tol
    warning('Frame rate is not constant...')
end

firstpulstime_lfp = pulset(1);

%% Reset time to align to LFP

t_align = t_EMG-firstpulstime_abf+firstpulstime_lfp;
Whints = [wh_on wh_off]-firstpulstime_abf+firstpulstime_lfp;
NWhints = [nwh_on nwh_off]-firstpulstime_abf+firstpulstime_lfp;

%% Figure
if SHOWFIG
    switch EMGanalysis
        case 'EMG'
      
            figure
            subplot(4,1,1)
            plot(t_align,EMGz,'k')

            hold on
            plot(t_align,EMGsm,'b','linewidth',2)
            plot(Whints',EMGparms.Whthreshold.*ones(size(Whints))','g','linewidth',2)
            plot(NWhints',EMGparms.NWhthreshold.*ones(size(NWhints))','r','linewidth',2)
            axis tight
            ylim([-100 100])
            ylabel('EMG (modZ)');

            subplot(4,1,2)
            plot(t_align,EMGz,'k')

            hold on
            plot(t_align,EMGsm,'b','linewidth',2)
            plot(Whints',EMGparms.Whthreshold.*ones(size(Whints))','g','linewidth',2)
            plot(NWhints',EMGparms.NWhthreshold.*ones(size(NWhints))','r','linewidth',2)
            xlim([100 160])
            ylim([-20 40])
            ylabel('EMG (modZ)');

            subplot(4,2,6)
            hist(log10(EMGsm),100)
            hold on
            plot([1 1].*log10(EMGparms.Whthreshold),get(gca,'ylim'),'g')
            plot([1 1].*log10(EMGparms.NWhthreshold),get(gca,'ylim'),'r')

            axis tight
            xlabel('EMG Envelope (modZ)');
            LogScale('x',10)
            xlim([-1.5 max(log10(EMGsm))])

            subplot(4,2,8)
            plot(durhist.bins,durhist.NWh,'r','linewidth',2)
            hold on
            plot(durhist.bins,durhist.Wh,'g','linewidth',2)
            LogScale('x',10)
            xlabel('Duration (s)')
            ylabel('# Epochs')
            legend('NWh','Wh')

            subplot(4,2,5)
            plot(t_abf-firstpulstime_abf+firstpulstime_lfp,pulse_abf,'k')
            hold on
            plot(pulset,pulsethreshold_abf.*ones(size(pulset)),'r+')
            xlim(firstpulstime_lfp+[-0.2 0.5])
            ylabel('Clampex Pulse Onset')

            subplot(4,2,7)

            plot(t_pulse,timepulses,'k')
            hold on
            plot(pulset,pulsethreshold.*ones(size(pulset)),'r+')
            xlim(firstpulstime_lfp+[-0.2 0.5])
            ylabel('Intan Pulse Onset')

            NiceSave('WhiskingDetection',figfolder,baseName)
        
        case 'touch'
            figure
            subplot(4,1,1)
            plot(t_align,EMGz,'k')

            hold on
            plot(t_align,EMGsm,'b','linewidth',2)
            plot(Whints',EMGparms.Whthreshold.*ones(size(Whints))','g','linewidth',2)
            plot(NWhints',EMGparms.NWhthreshold.*ones(size(NWhints))','r','linewidth',2)
            axis tight
            ylim([-100 100])
            ylabel('Piezo (modZ)');

            subplot(4,1,2)
            plot(t_align,EMGz,'k')

            hold on
            plot(t_align,EMGsm,'b','linewidth',2)
            plot(Whints',EMGparms.Whthreshold.*ones(size(Whints))','g','linewidth',2)
            plot(NWhints',EMGparms.NWhthreshold.*ones(size(NWhints))','r','linewidth',2)
            xlim([100 160])
            ylim([-20 40])
            ylabel('Piezo (modZ)');

            subplot(4,2,6)
            hist(log10(EMGsm),100)
            hold on
            plot([1 1].*log10(EMGparms.Whthreshold),get(gca,'ylim'),'g')
            plot([1 1].*log10(EMGparms.NWhthreshold),get(gca,'ylim'),'r')

            axis tight
            xlabel('Piezo Envelope (modZ)');
            LogScale('x',10)
            xlim([-1.5 max(log10(EMGsm))])

            subplot(4,2,8)
            plot(durhist.bins,durhist.NWh,'r','linewidth',2)
            hold on
            plot(durhist.bins,durhist.Wh,'g','linewidth',2)
            LogScale('x',10)
            xlabel('Duration (s)')
            ylabel('# Epochs')
            legend('NonTouch','Touch')

            subplot(4,2,5)
            plot(t_abf-firstpulstime_abf+firstpulstime_lfp,pulse_abf,'k')
            hold on
            plot(pulset,pulsethreshold_abf.*ones(size(pulset)),'r+')
            xlim(firstpulstime_lfp+[-0.2 0.5])
            ylabel('Clampex Pulse Onset')

            subplot(4,2,7)

            plot(t_pulse,timepulses,'k')
            hold on
            plot(pulset,pulsethreshold.*ones(size(pulset)),'r+')
            xlim(firstpulstime_lfp+[-0.2 0.5])
            ylabel('Intan Pulse Onset')

            NiceSave('TouchDetection',figfolder,baseName)
    end
end
%%

switch EMGanalysis
    case 'EMG'
        EMGwhisk.timestamps = t_align;
        EMGwhisk.EMG = EMGz;
        EMGwhisk.EMGsm = EMGsm;
        EMGwhisk.ints.Wh = Whints;
        EMGwhisk.ints.NWh = NWhints;
        EMGwhisk.samplingRate = sf_down;
        EMGwhisk.detectorparms = EMGparms;
        EMGwhisk.detectorname = 'GetWhiskFromEMG';
        EMGwhisk.detectiondate = today('datetime');
    case 'touch'
        Piezotouch.timestamps = t_align;
        Piezotouch.Piezo = EMGz;
        Piezotouch.Piezosm = EMGsm;
        Piezotouch.ints.Touch = Whints;
        Piezotouch.ints.NoTouch = NWhints;
        Piezotouch.samplingRate = sf_down;
        Piezotouch.detectorparms = EMGparms;
        Piezotouch.detectorname = 'GetWhiskFromEMG';
        Piezotouch.detectiondate = today('datetime');
end

if saveMat 
    switch EMGanalysis
        case 'EMG'
            save(savefile,'EMGwhisk')
        case 'touch'
            save(savefile,'Piezotouch')
    end
end
end

