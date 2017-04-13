function [ EMGwhisk ] = GetWhiskFromEMG( baseName,basePath )
%[ EMGwhisk ] = GetWhiskFromEMG( baseName,basePath ) 
%This is a detector that extracts whisking/nonwhisking epochs from 
%implanted EMG in the whisker pad. Extracts also the EMG and EMG envelope.
%
%INPUT
%Assumes presence of the following files:
%	basePath/baseName/baseName.abf    (whisking/camera pulses from clampex)
%   basePath/baseName/analogin.dat    (camera pulses to intan)
%If baseName,basePath are not specified as inputs, tries the current path.
%
%OUTPUT
%Creates file:
%   basePath/baseName/baseName.EMGwhisk.states.mat
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

if nargin==0
    [basePath,baseName] = fileparts(pwd);
end
%%
abfname = fullfile(basePath,baseName,[baseName,'.abf']);
analogName = fullfile(basePath,baseName,['analogin.dat']);
savefile = fullfile(basePath,baseName,[baseName,'.EMGwhisk.states.mat']);
figfolder = fullfile(basePath,baseName,'DetectionFigures');

if ~exist(abfname,'file')
    display('No .abf file named basePath/baseName/baseName.abf')
end
if ~exist(analogName,'file')
    display('No analogin.dat in basePath/baseName/')
end


%% Clampex File
timechan = 1;
emgchan = 2;
sf_abf = 20000; %Sampling Frequency of the .abf file


abffile = abfload(abfname);
pulse_abf = abffile(:,1);
EMG = abffile(:,2);

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
%%

EMGparms.gausswidth = 0.05;  %Gaussian width for smoothing (s)
EMGparms.Whthreshold = 3;    %EMG Threshold for Whisking (modSTDs)
EMGparms.NWhthreshold = 0.8;    %EMG Threshold for Whisking (modSTDs)
%EMGparms.threshold = 1;    %EMG Threshold for Whisking (modSTDs)
EMGparms.minwhisk = 0.1;     %Minimum whisking duration (s)
EMGparms.minNWh = 0.1;       %Minimum nonwhisking duration (s)
EMGparms.whiskmerge = 0.1;     %Minimum interwhisking duration (s)
EMGparms.NWhmerge = 0.02;       %Minimum internonwhisking duration (s)

%% Z-Score the EMGZ and get EMG envelope with RMS
EMGz = NormToInt(EMG,[],[],'modZ'); %Modified Z score - robust to outliers
EMGsm = RMSEnvelope(EMGz,EMGparms.gausswidth,1/sf_down);
EMGsm = EMGsm-min(EMGsm);

whisk.EMGenvelope = EMGsm;

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

%Merge brief interruptions
[ NWhints ] = MergeSeparatedInts( [nwh_on,nwh_off],EMGparms.NWhmerge );
[ Whints ] = MergeSeparatedInts( [wh_on,wh_off],EMGparms.NWhmerge );


%Drop nonwhisk epochs smaller than a minimum
[nwh_on,nwh_off] = MinEpochLength(NWhints(:,1),NWhints(:,2),EMGparms.minNWh,1/sf_down);
%Drop whisking epochs smaller than a minimum
[wh_on,wh_off] = MinEpochLength(Whints(:,1),Whints(:,2),EMGparms.minwhisk,1/sf_down);


%Convert to seconds
wh_on = wh_on./sf_down;
wh_off = wh_off./sf_down;
nwh_on = nwh_on./sf_down;
nwh_off = nwh_off./sf_down;

%% Durations
Whdur = wh_off-wh_on;
NWhdur = nwh_off-nwh_on;

numbinst = 25;
durhist.bins = linspace(-1,2,numbinst);
durhist.Wh = hist(log10(Whdur),durhist.bins);
durhist.NWh = hist(log10(NWhdur),durhist.bins);

%% Identify Pulses from Camera (in clampex)
pulsethreshold_abf =0.5;  %Adjust this later to set based on input.
pulseonsets = find(diff(pulse_abf<pulsethreshold_abf)==1);
pulset = t_abf(pulseonsets);

interpulse = diff(pulset);

%remove the first trigger, which is just camera onset... 
%Make this more rigorous later
if interpulse(1) > sum(interpulse(2:3))
    pulset(1) = []; 
end
firstpulstime_abf = pulset(1);

%% Load the analogin for the timestamps (pulses in intan)

timepulses = readmulti(analogName,1);

sf_pulse = 1./20000; %Sampling Frequency of the .abf file
t_pulse = [1:length(timepulses)]'.*sf_pulse;

pulsethreshold =1e4;  %Adjust this later to set based on input.
pulseonsets = find(diff(timepulses<pulsethreshold)==1);
pulset = t_pulse(pulseonsets);

minpulsedur = 0.003; %Remove double/noise crossings

shortpulses=diff(pulset)<(minpulsedur);
pulset(shortpulses) = [];

interpulse = diff(pulset);

if interpulse(1) > sum(interpulse(2:3))
    pulset(1) = []; 
end
firstpulstime_lfp = pulset(1);

%% Reset time to align to LFP

t_align = t_EMG-firstpulstime_abf+firstpulstime_lfp;
Whints = [wh_on wh_off]-firstpulstime_abf+firstpulstime_lfp;
NWhints = [nwh_on nwh_off]-firstpulstime_abf+firstpulstime_lfp;

%% Figure

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

%%

EMGwhisk.ints.Wh = Whints;
EMGwhisk.ints.NWh = NWhints;
EMGwhisk.detectorparms = EMGparms;
EMGwhisk.detectorname = 'GetWhiskFromEMG';
EMGwhisk.detectiondate = today;
EMGwhisk.EMG = EMGz;
EMGwhisk.EMGsm = EMGsm;
EMGwhisk.t = t_align;

save(savefile,'EMGwhisk')
end

