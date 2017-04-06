function [EMGCorr,sf_EMG] = EMGCorrForSleepscore(basenamepath,scoretime,specialchannels,specialshanks)
% Based on Erik Schomburg's work and code.  Grabs channels and calculates
% their correlations in the 300-600Hz band over sliding windows of 0.5sec.
% Channels are automatically selected and are a combination first channels
% on each shank (not using directly adjacent shanks, they must be at least 
% 2 shanks away.  
% First channels on "specialshanks" are always used
% "specialchannels" are also always used. 
% Special channels should be 0-indexed, per neuroscope convention
% Requires .lfp/lfp and .xml.  Assumes each spikegroup in the .xml
% represents a "shank"
% 
% Mean pairwise correlations are calculated for each time point.
% 
% Adapted 2015 Brendon Watson

%% Parameters
savebool = 0;

%% get basics about.lfp/lfp file
if strcmp(basenamepath(end-3:end),'.lfp') || strcmp(basenamepath(end-3:end),'.eeg')
    lfploc = basenamepath;
    xmlloc = [basenamepath(1:end-4),'.xml'];
    saveloc = [basenamepath(1:end-4),'_EMGCorr.mat'];
else
    if ~isempty(dir('*.lfp'))
       lfploc = [basenamepath '.lfp'];
    elseif ~isempty(dir('*.eeg'))
       lfploc = [basenamepath '.eeg'];
    else
        return
    end
    xmlloc = [basenamepath,'.xml'];
    saveloc = [basenamepath,'_EMGCorr.mat'];
end

Par = LoadParameters(xmlloc); % now using the updated version

Fs = Par.lfpSampleRate; % Hz, LFP sampling rate
nChannels = Par.nChannels;

if isfield(Par,'SpkGrps')
    SpkGrps = Par.SpkGrps;
elseif isfield(Par,'AnatGrps')
    SpkGrps = Par.AnatGrps;
    display('No SpikeGroups, Using AnatomyGroups')
elseif isfield(Par,'ElecGp')
    for i=1:length(Par.ElecGp.channel)
    SpkGrps.Channels = str2num(Par.ElecGp.channel);
else
    error('No SpikeGroups...')
end
    

xcorr_halfwindow_s = 0.5;%specified in s
% downsampleFs = 125;
% downsampleFactor = round(Fs/downsampleFs);
binScootS = 0.5;
sf_EMG = 1/binScootS;
binScootSamps = Fs*binScootS;
corrChunkSz = 20;%for batch-processed correlations


%% input handling: channel selection
if ~exist('specialchannels','var')
    specialchannels = [];
end
if ~exist('specialshanks','var')
    specialshanks = [];
end

[PATHSTR] = fileparts(basenamepath);

rejectchannels = [];
if exist(fullfile(PATHSTR,'bad_channels.txt'),'file')%bad channels is an ascii/text file where all lines below the last blank line are assumed to each have a single entry of a number of a bad channel (base 0)
    t = ReadBadChannels_ss(PATHSTR);
    rejectchannels = cat(1,rejectchannels(:),t(:));
end



%% Pick shanks to analyze
% get spike groups,
% pick every other one... unless specialshanks, in which case pick non-adjacent
%This is potentially dangerous in combination with rejectchannels... i.e.
%what if you pick every other shank but then the ones you pick are all
%reject because noisy shank.
for nn = 1:2  %This loop is part of the crappy fix for that problem

spkgrpstouse = nn:2:length(SpkGrps);
spkgrpstouse = spkgrpstouse(:);

if ~isempty(specialshanks) 
    spkgrpstouse = unique(cat(1,spkgrpstouse,specialshanks));
end
if ~isempty(specialchannels)
    for a = 1:length(specialchannels)
        for b = 1:length(SpkGrps)
            if ismember(specialchannels(a),SpkGrps(b).Channels)
                spkgrpstouse = cat(1,spkgrpstouse,b);
                continue
            end
        end
    end
end

%Go through each spike group and add one channel to the list of channels to
%use for calculating high-f correlation.
xcorr_chs = [];
for a = 1:length(spkgrpstouse)
    spkgroup = spkgrpstouse(a);
    %Are any of the channels in this spike group "special channels"?
    [lia,lib]=ismember(specialchannels,SpkGrps(spkgroup).Channels);
    if any(lia)
        xcorr_chs(end+1) = SpkGrps(spkgroup).Channels(lib);
    else
        availchans = SpkGrps(spkgroup).Channels;
        %Don't add rejectchannels to the list of available channels
        availchans = setdiff(availchans,rejectchannels); 

        
        if length(availchans)>=1
            xcorr_chs(end+1) = availchans(1); %Pick the first channel
        end
    end
end

xcorr_chs = unique(xcorr_chs);

%This is a crappy catch for the issue of no channels due to too many
%rejectchannels
if length(xcorr_chs)>=2
    break
elseif length(xcorr_chs)<2 && nn==2
    display('You have no channels for EMG... This is a bug?')
    return
end
end
    

%% Read and filter channel
% read channels
xcorr_chs = xcorr_chs + 1; % loadparameters returns 0 indexed (neuroscope) channels, 
                           % but Loadbinary.m takes 1-indexed channel #'s
lfp = LoadBinary(lfploc,'nChannels',nChannels,'channels',xcorr_chs); %read and convert to mV    
% Filter first in high frequency band to remove low-freq physiologically
% correlated LFPs (e.g., theta, delta, SPWs, etc.)

% eeg = LoadBinary_Down_ss(lfploc,'frequency',Fs,...
%     'nchannels',nChannels,'channels',xcorr_chs+1,...
%     'start',scoretime(1),'duration',diff(scoretime));
% %+1 is applied to channel numbers here for 0 (neuroscope) vs 1 (LoadBinary)
% %indexing.


xcorr_freqband = [275 300 600 625]; % Hz
lfp = filtsig_in(lfp, Fs, xcorr_freqband);

%% xcorr 'strength' is the summed correlation coefficients between channel
% pairs for a sliding window of 25 ms
xcorr_window_samps = round(xcorr_halfwindow_s*Fs);
xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;%+- that number of ms in samples

% xcorr 'strength' is the averaged correlation coefficients between channel
% old version... single correlation calculated at once
% % pairs for a sliding window of 25 ms
% xcorr_window_samps = round(xcorr_window_s*Fs);
% xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;%+- that number of ms in samples
% 
% numbins = (size.lfp,1) - xcorr_window_samps*2)/binScootSamps;
% xcorrStrength = zeros(numbins, 1);
% binind = 0;
% for i=(1+xcorr_window_inds(end)):binScootSamps:(size.lfp,1)-xcorr_window_inds(end))
%     binind = binind+1;
%     for j=1:(length(xcorr_chs)-1)
%         for k=(j+1):length(xcorr_chs)
%             s1 =.lfp(i + xcorr_window_inds, j);
%             s2 =.lfp(i + xcorr_window_inds, k);
%             tmp = corrcoef(s1,s2);
%             xcorrStrength(binind) = xcorrStrength(binind) + tmp(1,2);
%         end
%     end
% end
% xcorrStrength = xcorrStrength/(length(xcorr_chs)*(length(xcorr_chs)-1)/2);


% new version... batches of correlation calculated at once
timestamps = (1+xcorr_window_inds(end)):binScootSamps:(size(lfp,1)-xcorr_window_inds(end));
numbins = length(timestamps);
EMGCorr = zeros(numbins, 1);
% tic
counter = 1;
for j=1:(length(xcorr_chs)-1)
    for k=(j+1):length(xcorr_chs)
        disp(counter)
        c1 = [];
        c2 = [];
        binind = 0;
        binindstart = 1;
        for i = timestamps
            binind = binind+1;
            s1 =lfp(i + xcorr_window_inds, j);
            s2 =lfp(i + xcorr_window_inds, k);
            c1 = cat(2,c1,s1);
            c2 = cat(2,c2,s2);
            if size(c1,2) == corrChunkSz || i == timestamps(end)
                binindend = binind;
                tmp = corr(c1,c2);
                tmp = diag(tmp);
                EMGCorr(binindstart:binindend) = EMGCorr(binindstart:binindend) + tmp;
                c1 = [];
                c2 = [];
                binindstart = binind+1;
            end
            counter = counter+1;
        end
    end
end
% toc

EMGCorr = EMGCorr/(length(xcorr_chs)*(length(xcorr_chs)-1)/2);

EMGCorr = cat(2,timestamps'/Fs,EMGCorr);
ChannelsCompared = xcorr_chs;
% EMGCorrData = v2struct_ss(EMGCorr,ChannelsCompared,AnatShankSite);
EMGCorrData = v2struct(EMGCorr,ChannelsCompared);
if savebool
    % save...
    save(saveloc,'EMGCorrData','EMGCorr','sf_EMG');
end

function [filt_sig, Filt] = filtsig_in(sig, Fs, filtband_or_Filt)
% [filt_sig, Filt] = filtsig(sig, dt_ms, filtband_or_Filt)
%
% Created by: Erik Schomburg, 2011

if isnumeric(filtband_or_Filt)
    h  = fdesign.bandpass(filtband_or_Filt(1), filtband_or_Filt(2), filtband_or_Filt(3), filtband_or_Filt(4), ...
        60, 1, 60, Fs);
    Filt = design(h, 'butter', 'MatchExactly', 'passband');
else
    Filt = filtband_or_Filt;
end

if ~isempty(sig)
    if iscell(sig)
        filt_sig = cell(size(sig));
        for i=1:length(sig(:))
            filt_sig{i} = filter(Filt, sig{i});
            filt_sig{i} = filter(Filt, filt_sig{i}(end:-1:1));
            filt_sig{i} = filt_sig{i}(end:-1:1);
        end
    elseif ((size(sig,1) > 1) && (size(sig,2) > 1))
        filt_sig = zeros(size(sig));
        for i=1:size(filt_sig,2)
            filt_sig(:,i) = filter(Filt, sig(:,i));
            filt_sig(:,i) = filter(Filt, filt_sig(end:-1:1,i));
            filt_sig(:,i) = filt_sig(end:-1:1,i);
        end
    else
        filt_sig = filter(Filt, sig);
        filt_sig = filter(Filt, filt_sig(end:-1:1));
        filt_sig = filt_sig(end:-1:1);
    end
else
    filt_sig = [];
end

