function [ EMGwhisk ] = GetWhiskFromPiezo( basePath,varargin )
%[ EMGwhisk ] = GetWhiskFromEMG(basePath ) 
%This is a detector that extracts whisking/nonwhisking epochs from 
%implanted EMG in the whisker pad. Extracts also the EMG and EMG envelope.
%
%INPUT
%   Assumes presence of the following files:
%       basePath/baseName.abf    (whisking/camera pulses from clampex)
%       basePath/analogin.dat    (camera pulses to intan)
%   where basePath is a folder of the form: 
%       whateverPath/baseName/
%If basePath not specified, tries the current path.
%
%   (options)
%       'PulseChannel'
%       'EMGChannel'
%
%OUTPUT
%Creates file:
%   basePath/baseName.EMGwhisk.states.mat
%
%
%
%To Add: option for no camera pulses, which just align time to intan/pupil
%
%
%DLevenstein 2017
%% DEV
% basePath = '/mnt/proraidDL/Database/WMProbeData/';
% baseName = 'Layers_LFP_Test02_170323_151411';
%%

if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);
%%
abfname = fullfile(basePath,[baseName,'.abf']);
analogName = fullfile(basePath,['analogin.dat']);
savefile = fullfile(basePath,[baseName,'.EMGwhisk.states.mat']);
figfolder = fullfile(basePath,'DetectionFigures');

if ~exist(abfname,'file')
    display('No .abf file named basePath/baseName/baseName.abf')
end
if ~exist(analogName,'file')
    display('No analogin.dat in basePath/baseName/')
end


%% Clampex File

[abffile,si,file_info] = abfload(abfname);
sf_abf = 1./(si.*10^-6);

%Prompt user for channel
numabfchans = length(abffile(1,:));
if ~exist('chanNums','var')
    piezochan = listdlg('ListString',file_info.recChNames,...
        'PromptString',['Which channels is Whisker? ']);
    %Replace this with prompt file_info.recChNames
end
%%
EMG = abffile(:,piezochan);
t_abf = [1:length(EMG)]'./sf_abf;

%%
% figure
% plot(t_abf,EMG,'k')
%% Get the whisking envelope
%EMGrange = [400 3000];
%EMG = FiltNPhase(EMG,EMGrange,sf_abf);

downsamplefactor = 16; %Downsample to same as the LFP;
EMG = downsample(EMG,downsamplefactor);
t_EMG = downsample(t_abf,downsamplefactor);
sf_down = sf_abf./downsamplefactor;
%%

EMGparms.gausswidth = 0.05;  %Gaussian width for smoothing (s)
EMGparms.Whthreshold = 3;    %EMG Threshold for Whisking (modSTDs)
EMGparms.NWhthreshold = 0.5;    %EMG Threshold for Whisking (modSTDs)
%EMGparms.threshold = 1;    %EMG Threshold for Whisking (modSTDs)
EMGparms.minwhisk = 0.1;     %Minimum whisking duration (s)
EMGparms.minNWh = 0.1;       %Minimum nonwhisking duration (s)
EMGparms.whiskmerge = 0.1;     %Minimum interwhisking duration (s)
EMGparms.NWhmerge = 0.02;       %Minimum internonwhisking duration (s)

%% Z-Score the EMGZ and get EMG envelope with RMS
EMGz = NormToInt(EMG,'modZ'); %Modified Z score - robust to outliers
EMGsm = RMSEnvelope(EMGz,EMGparms.gausswidth,1/sf_down);
EMGsm = EMGsm-min(EMGsm);
EMGwhisk.EMGenvelope = EMGsm;

%% Set the thresholds by whisking troughs - 
%find by "gradient descent"(ish) from initial guess
EMGbins = linspace(-1.5,2,100);
EMGhist = hist(log10(EMGsm),EMGbins);
EMGgrad = smooth(gradient(EMGhist),4);

%Find troughs (gradient crossing from - to +)
troughidx = find(diff(EMGgrad>0)==1);
troughs = 10.^EMGbins(troughidx);

%Get sign of gradient at each of the thresholds and use that to pick trough
Whsign = sign(interp1(EMGbins,EMGgrad,log10(EMGparms.Whthreshold),'nearest'));
if Whsign==-1
    EMGparms.Whthreshold = troughs(find(troughs>EMGparms.Whthreshold,1,'first'));
    if isempty(EMGparms.Whthreshold)
        EMGparms.Whthreshold = troughs(end);
    end
elseif Whsign==1
    EMGparms.Whthreshold = troughs(find(troughs<EMGparms.Whthreshold,1,'last'));
end

NWhsign = sign(interp1(EMGbins,EMGgrad,log10(EMGparms.NWhthreshold),'nearest'));
if NWhsign==-1
    EMGparms.NWhthreshold = troughs(find(troughs>EMGparms.NWhthreshold,1,'first'));
elseif NWhsign==1
    EMGparms.NWhthreshold = troughs(find(troughs<EMGparms.NWhthreshold,1,'last'));
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



%% Reset time to align to LFP

%t_EMG = t_EMG-firstpulstime_abf+firstpulstime_lfp;
Whints = [wh_on wh_off];
NWhints = [nwh_on nwh_off];

%% Figure

figure
subplot(4,1,1)
    plot(t_EMG,EMGz,'k')

    hold on
    plot(t_EMG,EMGsm,'b','linewidth',2)
    plot(Whints',EMGparms.Whthreshold.*ones(size(Whints))','g','linewidth',2)
    plot(NWhints',EMGparms.NWhthreshold.*ones(size(NWhints))','r','linewidth',2)
    axis tight
    ylim([-100 100])
    ylabel('EMG (modZ)');
    
subplot(4,1,2)
    plot(t_EMG,EMGz,'k')

    hold on
    plot(t_EMG,EMGsm,'b','linewidth',2)
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


NiceSave('WhiskingDetection',figfolder,baseName)

%%

EMGwhisk.ints.Wh = Whints;
EMGwhisk.ints.NWh = NWhints;
EMGwhisk.detectorparms = EMGparms;
EMGwhisk.detectorname = 'GetWhiskFromEMG';
EMGwhisk.detectiondate = today('datetime');
EMGwhisk.EMG = EMGz;
EMGwhisk.EMGsm = EMGsm;
EMGwhisk.t = t_EMG;

save(savefile,'EMGwhisk')
end

