function [ pSpindleInts,cycletimemap,deltapeaks,SpindleStats ] = FindSpindlesAndSWs(datasetfolder,recname,figfolder,StateIntervals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%




%% LOAD THE SESSION IN FMA
%recname = 'c3po_160202';
%recname = '20140526_277um';
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/~updated/Recordings (1)/';
%recname = '20140526_277um';
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/DTData/';
%recname = 'DT2_rPPC_rCCG_362um_218um_20160209_160209_183610';
%figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisFigsMats/SpindleID/';
% probegroups = {1:5,7:12}; %Group 1: Cingulate.  %Group2: Parietal DT2, spike group 6,13 are top of shank
% %probegroups = {7:12,1:6}; %Group 1: Cingulate.  %Group2: HPC c3po
% probegroups = {1:4,5:6}; %Group 1: Cingulate.  %Group2: HPC jenn1
%numprobes = length(probegroups);



xmlfilename = [datasetfolder,'/',recname,'/',recname,'.xml'];
if exist (fullfile(datasetfolder,recname,[recname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recname,[recname,'.lfp']);
elseif exist (fullfile(datasetfolder,recname,[recname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recname,[recname,'.lfp']);
else 
    display('No .lfp file')
end


SetCurrentSession(xmlfilename)
global DATA;

spikegroups = DATA.spikeGroups.groups;

numsites = DATA.nChannels;
numgroups = length(spikegroups);


% Par = LoadPar(xmlfilename);
% Fs = Par.lfpSampleRate; % Hz, LFP sampling rate
% nChannels = Par.nChannels;
% spikegroups = {Par.SpkGrps(:).Channels};
numgroups = length(spikegroups);

%% Load the channel map to get the cortical channels
CTXlabels = {'mPFC','ACC','MotorCtx','OFC'};

spikegroupanatomyfilename = fullfile(datasetfolder,recname,[recname,'_SpikeGroupAnatomy.csv']);
spkgroupanatomy=readtable(spikegroupanatomyfilename);

ctxgroups = ismember(spkgroupanatomy.AnatomicalSite,CTXlabels);
ctxgroups = spkgroupanatomy.SpikeGroup(ctxgroups);
ctxchannels = [spikegroups{ctxgroups}];

%%

if ~exist('StateIntervals','var')
    load(fullfile(datasetfolder,recname,[recname,'_SleepScore.mat']));
end


NREMint = StateIntervals.NREMstate;

if exist (fullfile(datasetfolder,recname,[recname,'_GoodSleepInterval.mat']),'file')
   load(fullfile(datasetfolder,recname,[recname,'_GoodSleepInterval.mat']));
   [NREMint] = RestrictInts(NREMint,GoodSleepInterval); 
end


%% Identidy Spindles

[Spindles,ChannelProperties] = DetectSPINDLES(rawlfppath,ctxchannels,NREMints,'best',figfolder);


%% Adjust Start End Time to peaks, return spindle cycle normalized time
%Put this into DetectSPINDLES

%% Delta Peak Times
[deltapeaks] = DeltaPeakTimes(ctxchannels,NREMint,figfolder,recname);


%% Characteriation of pSpindle LFP and peaks (add peaks characterization from below)
[SpindleStats] = SpindleIntLFP(pSpindleInts,cycletimemap,deltapeaks,ctxchannels,NREMint,figfolder,recname)



end

