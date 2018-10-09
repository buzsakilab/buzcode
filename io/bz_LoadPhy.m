
function [spikes] = bz_LoadPhy(varargin)
% Load kilosort/phy clusters
%
% USAGE
%
%   [spikes] = bz_loadPhy(varargin)
%
% INPUTS
%
% basepath        -path to recording (where .dat are)
% kilosort_path   -path to kilosort* folder after phy curing (by default
%                   .\kilosort_\)
% getWaveforms    -logical (default=true) to load mean of raw waveform data
% saveMat         -logical (default=false) to save in buzcode format
% UID             -vector subset of UID's to load 
% fs              -scalar (default is taken from sessionInfo, otherwise
%                   30000). Sampling frequency. 
% nChannels       -scalar (default is taken from sessionInfo, otherwise
%                   32). Total number of channels.
% forceReload    -logical (default=false) to force loading from
%                   res/clu/spk files
% getFeatures    -logical (default=true) to compute cluster features
% showFeatures   -logical (default=false) to show cluster features for each
%                   cluster. This option force getFeatures to be true.
%
% OUTPUTS
%
% spikes - cellinfo struct with the following fields
%   .sessionName    -name of recording file
%   .UID            -unique identifier for each neuron in a recording
%   .times          -cell array of timestamps (seconds) for each neuron
%   .spindices      -sorted vector of [spiketime UID], useful for 
%                           input to some functions and plotting rasters
%   .shankID        -shank ID that each neuron was recorded on
%   .amplitude:     amplitudes for cluster N
%   .maxWaveformCh  -channel # with largest amplitude spike for each neuron
%   .rawWaveform    -average waveform on maxWaveformCh (from raw .dat)
%
%   MV 2018

% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'kilosort_path',ls('Kilosort*'),@isstr); % probably this line only works in windows
addParameter(p,'getWaveforms',true,@islogical)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'UID',[],@isvector);
addParameter(p,'fs',30000,@isnumeric);
addParameter(p,'nChannels',32,@isnumeric);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'getFeatures',true,@islogical);
addParameter(p,'showFeatures',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
kilosort_path = p.Results.kilosort_path;
getWave = p.Results.getWaveforms;
saveMat = p.Results.saveMat;
UID = p.Results.UID;
fs = p.Results.fs; % it will be overwritten if bz_getSessionInfo
nChannels = p.Results.nChannels; % it will be overwritten if bz_getSessionInfo
forceReload = p.Results.forceReload;
getFeat = p.Results.getFeatures;
showFeat = p.Results.showFeatures;

if showFeat; getFeat = true; end

try [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', false);
    fs = sessionInfo.rates.wideband;
    nChannels = sessionInfo.nChannels;
catch
    warning('SessionInfo file not found.');
end

if exist([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'],'file') && forceReload == false
    disp('loading spikes from cellinfo file..')
    load([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'])
    
else
    spike_cluster_index = readNPY(fullfile(kilosort_path, 'spike_clusters.npy'));
    spike_times = readNPY(fullfile(kilosort_path, 'spike_times.npy'));
    % spike_amplitudes = readNPY(fullfile(kilosort_path, 'amplitudes.npy'));
    % spike_clusters = unique(spike_cluster_index);
    cluster_group = tdfread(fullfile(kilosort_path,'cluster_group.tsv'));
    shanks = readNPY(fullfile(kilosort_path, 'shanks.npy')); % done


    spikes = [];
    spikes.sessionName = sessionInfo.FileName;
    jj = 1;
    for ii = 1:length(cluster_group.group)
        if strcmp(cluster_group.group(ii,:),'good ')
            ids = find(spike_cluster_index == cluster_group.cluster_id(ii)); % cluster id
            spikes.UID(jj) = cluster_group.cluster_id(ii);
            spikes.times{jj} = double(spike_times(ids))/fs; % cluster time
            spikes.ts{jj} = double(spike_times(ids)); % cluster time
            cluster_id = find(cluster_group.cluster_id == spikes.UID(jj));
            spikes.shankID(jj) = shanks(cluster_id);
            % spikes.amplitudes{jj} = double(spike_amplitudes(ids));
            jj = jj + 1;
        end
    end

    % get waveforms
    nPull = 500; % number of spikes to pull out
    wfWin = 0.003; % Size of waveform windows
    hpFilt = designfilt('highpassiir','FilterOrder',8, 'PassbandFrequency',1000,'PassbandRipple',0.1, 'SampleRate',fs);

    f = waitbar(0,'Getting waveforms...');
    wfWin = round((wfWin * fs)/2);
    if getWave
        for ii = 1 : size(spikes.times,2)
            spkTmp = spikes.times{ii};
            if length(spkTmp) > nPull
                spkTmp = spkTmp(randperm(length(spkTmp)));
                spkTmp = spkTmp(1:nPull);
            end
            wf = [];
            for jj = 1 : length(spkTmp)
                wf = cat(3,wf,bz_LoadBinary([sessionInfo.session.name '.dat'],'offset',spikes.ts{ii}(jj) - wfWin,...
                    'samples',wfWin * 2,'frequency',sessionInfo.rates.wideband,'nChannels',sessionInfo.nChannels));
            end
            wf = mean(wf,3);
            for jj = 1 : size(wf,2)
                wfF(:,jj) = filtfilt(hpFilt,wf(:,jj));
            end
            [~, spikes.maxWaveformCh(ii)] = max(abs(wfF(wfWin,:)));
            spikes.rawWaveform{ii} = wf(:,spikes.maxWaveformCh(ii)) - mean(wf(:,spikes.maxWaveformCh(ii)));
%            spikes.allWf{ii} = wfF;
            
            waitbar(ii/size(spikes.times,2),f,'Pulling out waveforms...');
        end
        close(f)
    end
    
    % spike measures
    disp('Computing spike features... ');
    if showFeat; figure; end
    if getFeat
        for ii = 1 : size(spikes.times,2)
            % firing rate
            spikes.firing_rate(ii) = size(spikes.times{ii},1)/(spikes.times{ii}(end) - spikes.times{ii}(1));
            
            % peak (neg) to peak (pos) duration, as Senzai et al 2017
            [~,tmp] = max(spikes.rawWaveform{ii}(size(spikes.rawWaveform{ii},1)/2:end));
            spikes.spk_duration(ii) = (tmp)/fs; % peak (negative) to peak (second positive) duration
            tmp = tmp + size(spikes.rawWaveform{ii},1)/2 - 1;
            
            % half width
            interpFac = 50;
            mean_spike = interp1(linspace(-wfWin,wfWin,length(spikes.rawWaveform{ii})), ...
                spikes.rawWaveform{ii} - mean(spikes.rawWaveform{ii}),...
                linspace(-wfWin,wfWin,interpFac * length(spikes.rawWaveform{ii})));
            [~,cutpoint_1] = min(abs(mean_spike(1 : wfWin * interpFac) - (mean_spike(wfWin * interpFac)/2)));
            [~,cutpoint_2] = min(abs(mean_spike(wfWin * interpFac + 1 : end) - (mean_spike(wfWin * interpFac)/2)));
            cutpoint_2 = cutpoint_2 + wfWin * interpFac;
            spikes.half_width(ii) = ((cutpoint_2 - cutpoint_1)/ interpFac)/fs;
            
            % asymmetry
            spkTemp = spikes.rawWaveform{ii};
            [pks, locs] = findpeaks(spkTemp);
            if isempty(locs(locs < size(spikes.rawWaveform{ii},1)/2)) 
                [~,peak1Loc] = min(abs(spkTemp(1:size(spkTemp,1)/2))); 
            else
                peak1Loc = locs(find(max(locs(locs < size(spkTemp,1)/2))==locs));
            end
            peak1 = spkTemp(peak1Loc);
            
            if isempty(pks(locs > size(spkTemp,1)/2))
                peak2Loc = size(spkTemp,1);
            else
                peak2Loc = locs(find(max(pks(locs > size(spkTemp,1)/2)) == pks));
            end
            peak2 = spkTemp(peak2Loc);
                
            spikes.asymmetry(ii) = (peak2 - peak1)/ (peak1 + peak2);
            
            % AUTOCORRELOGRAM FEATURES & DOUBLE EXPONENTIAL FITTING MODEL
            try ACG_mat = 1000*(CrossCorr(spikes.times{ii},spikes.times{ii},.001,100)/length(spikes.times{ii}));
                ACG_mat(51) = 0;
                [fmodel,~,~,paut] = fitpyrint(ACG_mat',0:50,0,20);
                spikes.ACG.fmodel{ii} = fmodel;
                spikes.ACG.ydata{ii} = ACG_mat;
                spikes.ACG.xdata{ii} = linspace(-0.050,0.05, length(ACG_mat));
                spikes.doubleExponentialACG(ii,:) =  paut;
            catch 
                warning('CrossCorr and fitpyrint not found. ACG can not be computed! ');
            end
            
            % Burstiness (As in Mizuseki et al, 2011). Fraction of spikes
            % with a ISI for following or preceding spikes < 0.006
            
            bursty = [];
            for jj = 2 : length(spikes.times{ii}) - 1
                bursty(jj) =  any(diff(spikes.times{ii}(jj-1 : jj + 1)) < 0.006);
            end 
            spikes.burstIndex(ii) = length(find(bursty > 0))/length(bursty);
            
            % plot features
            if showFeat
                xax = linspace(-wfWin/fs * 1000, wfWin/fs * 1000, length(spikes.rawWaveform{ii}));
                subplot(1,2,1)
                cla
                hold on
                plot(xax, spikes.rawWaveform{ii}); axis tight
                plot([xax(peak1Loc) xax(peak1Loc)], [0 spikes.rawWaveform{ii}(peak1Loc)],'LineWidth',1.5);
                plot([xax(peak2Loc) xax(peak2Loc)], [0 spikes.rawWaveform{ii}(peak2Loc)],'LineWidth',1.5);
                
                plot([0 xax(tmp)],[spikes.rawWaveform{ii}(tmp) spikes.rawWaveform{ii}(tmp)],'LineWidth',1.5);
                plot([xax(size(spikes.rawWaveform{ii},1)/2)],...
                    spikes.rawWaveform{ii}(size(spikes.rawWaveform{ii},1)/2),'o', 'LineWidth',1.5);
                xlabel('ms'); ylabel('amp');
                try subplot(1,2,2)
                    cla
                    area(spikes.ACG.xdata{ii}, spikes.ACG.ydata{ii},'LineStyle','none');
                    xlabel('ms'); ylabel('#');
                end
                
                disp('----------');
                fprintf('Freq rate: %6.4f \n', spikes.firing_rate(ii));
                fprintf('Asym: %6.4f \n', spikes.asymmetry(ii));
                fprintf('Burstiness: %6.4f \n', spikes.burstIndex(ii));
                fprintf('Spk duration: %6.4f \n\n', spikes.spk_duration(ii));
                spikes.label(ii) = input('Any key to continue... ');
                disp('----------');
                
            end
        end
    end

end

% saveMat (only saving if no exclusions)
if saveMat
    save([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'],'spikes');
end

% filter by UID input
if ~isempty(UID)
      [toRemove] = ~ismember(spikes.UID,UID);
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
         % spikes.region{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    %spikes.region = removeEmptyCells(spikes.region);
    % spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    
    if getWaveforms
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.rawWaveform{r} = [];
        end
    end
    spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
    spikes.maxWaveformCh(toRemove) = [];
    end
end

% Generate spindices matrics
for cc = 1:length(spikes.UID)
    groups{cc}=spikes.UID(cc).*ones(size(spikes.times{cc}));
end
if ~isempty(spikes.UID)
    alltimes = cat(1,spikes.times{:}); groups = cat(1,groups{:}); %from cell to array
    [alltimes,sortidx] = sort(alltimes); groups = groups(sortidx); %sort both
    spikes.spindices = [alltimes groups];
end

end