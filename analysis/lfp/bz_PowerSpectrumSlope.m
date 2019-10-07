function [specslope,spec] = bz_PowerSpectrumSlope(lfp,winsize,dt,varargin)
%[specslope,spec] = bz_PowerSpectrumSlope(lfp,winsize,dt) calculates the
%slope of the power spectrum, a metric of cortical state and E/I balance
%see Gao, Peterson, Voytek 2016;  Waston, Ding, Buzsaki 2017
%
%INPUTS
%   lfp         a buzcode-formatted lfp structure (use bz_GetLFP)
%               needs fields: lfp.data, lfp.timestamps, lfp.samplingRate
%   winsize     size of the silding time window (s, 2-4 recommended)
%   dt          sliding time interval (s)
%
%   (optional)
%       'frange'    (default: [4 100])
%       'channels'  subset of channels to calculate PowerSpectrumSlope
%                   (default: all)
%       'ints'      Intervals in which to calculate PSS. Default:[0 Inf]
%       'showfig'   true/false - show a summary figure of the results
%                   (default:false)
%       'saveMat'   put your basePath here to save/load
%                   baseName.PowerSpectrumSlope.lfp.mat  (default: false)
%       'Redetect'  (default: false) to force redetection even if saved
%                   file exists
%       'IRASA'     (default: true) use IRASA method to median-smooth power
%                   spectrum before fitting
%                   (Muthukumaraswamy and Liley, NeuroImage 2018)
%       'nfreqs'    number of frequency values used to fitting (default: 200)
%       'spectype'  'fft' or 'wavelet' (default: fft)
%                   if using wavelets, winsize corresponds to number of
%                   cycles (recommended: 10) and dt downsamples wavelets to
%                   approximate desired dt (recommended <0.01 for
%                   resolution up to 100Hz)
%                   
%
%OUTPUTS
%   specslope
%       .data           [Nt x NChannels] vector of the power spectrum slope   
%       .timestamps     [Nt] timestamps
%       .intercept      [Nt x NChannels] vector of the 0-intercept
%       .specgram       [Nt x Nfreqs x Nchannels] the spectrogram
%       .resid          [Nt x Nfreqs x Nchannels] the PSS-removed spectrogram
%   specgram
%       .data       complex-valued spectrogram
%       .timestamps
%       .amp        log10-transformed amplitude of the spectrogram
%
%
%DLevenstein 2018
%with IRASA code modified from R. Hardstone & W. Munoz - 2019
%%
p = inputParser;
addParameter(p,'showfig',false,@islogical)
addParameter(p,'saveMat',false)
addParameter(p,'channels',[])
addParameter(p,'frange',[4 100])
addParameter(p,'Redetect',false)
addParameter(p,'IRASA',true)
addParameter(p,'nfreqs',200)
addParameter(p,'spectype','fft')
addParameter(p,'ints',[0 inf])
parse(p,varargin{:})
SHOWFIG = p.Results.showfig;
saveMat = p.Results.saveMat;
channels = p.Results.channels;
frange = p.Results.frange;
REDETECT = p.Results.Redetect;
IRASA = p.Results.IRASA;
nfreqs = p.Results.nfreqs;
spectype = p.Results.spectype;
ints = p.Results.ints;

%%
if saveMat
    basePath = saveMat;
    baseName = bz_BasenameFromBasepath(basePath);
    savename = fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']);
    
    if exist(savename,'file') && ~REDETECT
        load(savename)
        return
    end
end

%% For multiple lfp channels
if ~isempty(channels)
    usechans = ismember(lfp.channels,channels);
    lfp.data = lfp.data(:,usechans);
    lfp.channels = lfp.channels(usechans);
elseif ~isfield(lfp,'channels')
    lfp.channels = nan;
end

if length(lfp.channels)>1
  %loop each channel and put the stuff in the right place
	for cc = 1:length(lfp.channels)
        bz_Counter(cc,length(lfp.channels),'Channels');
        lfp_temp = lfp; 
        lfp_temp.data = lfp_temp.data(:,cc); 
        lfp_temp.channels = lfp_temp.channels(cc);
        specslope_temp = bz_PowerSpectrumSlope(lfp_temp,winsize,dt,varargin{:},'saveMat',false);
        
        if ~exist('specslope','var')
            specslope = specslope_temp;
            %Don't save the big stuff for multiple channels.
            %specslope.resid = [];
        else
            specslope.data(:,cc) = specslope_temp.data;
            specslope.intercept(:,cc) = specslope_temp.intercept;
            specslope.resid(:,:,cc) = specslope_temp.resid;
            specslope.rsq(:,:,cc) = specslope_temp.rsq;
            specslope.specgram(:,:,cc) = specslope_temp.specgram;
        end
        
	end
    specslope.channels = lfp.channels;
    %spec = [];
    if saveMat
        save(savename,'specslope','-v7.3');
    end
    return
end
%% Calcluate spectrogram
noverlap = winsize-dt;

if IRASA
    maxRescaleFactor = 2.9; %as per Muthukumaraswamy and Liley, NeuroImage 2018
    %Add the padding for edge effects of frequency smoothings (IRASA)
    padding = floor((nfreqs./2).*log10(maxRescaleFactor^2)./log10(frange(2)./frange(1)));
    actualRescaleFactor = 10.^((log10(frange(2)./frange(1)).*padding)./(nfreqs-1)); %Account for rounding
    nfreqs = nfreqs+2*padding;
    frange = frange .* [1/actualRescaleFactor actualRescaleFactor];
end

switch spectype
    case 'fft'
        spec.freqs = logspace(log10(frange(1)),log10(frange(2)),nfreqs);
        winsize_sf = round(winsize .*lfp.samplingRate);
        noverlap_sf = round(noverlap.*lfp.samplingRate);
        [spec.data,~,spec.timestamps] = spectrogram(single(lfp.data),...
            winsize_sf,noverlap_sf,spec.freqs,lfp.samplingRate);
        
        %Interpolate the time stamps to match the LFP timestamps
        assumedLFPtimestamps = [0:length(lfp.data)-1]./lfp.samplingRate;
        spec.timestamps = interp1(assumedLFPtimestamps,lfp.timestamps,spec.timestamps,'nearest');
        
        keeptimes = InIntervals(spec.timestamps,ints);
        
        spec.amp = log10(abs(spec.data(:,keeptimes)'));
        spec.data = spec.data(:,keeptimes)';
        spec.timestamps = spec.timestamps(keeptimes)';
    case 'wavelet'
        downsampleout = round(lfp.samplingRate.*dt);
        dt = downsampleout/lfp.samplingRate; %actual dt another option is to do movmean on the amplitude...
        [spec] = bz_WaveSpec(lfp,'frange',frange,'nfreqs',nfreqs,'ncyc',winsize,...
            'downsampleout',downsampleout,'intervals',ints);
        
        spec.amp = log10(abs(spec.data));
end



%% IRASA before fitting
if IRASA
    %Frequencies inside padding 
    validFreqInds = padding+1 : nfreqs-padding;
    
    %Median smooth the spectrum
    for i_freq = validFreqInds
        inds = [i_freq-padding:i_freq-1 i_freq+1:i_freq+padding];
        resampledData(:,i_freq) = nanmedian(spec.amp(:,inds),2);
    end
    
    %Calculate the residual from the smoothed spectrum
    %spec.osci = (spec.amp(:,validFreqInds))-(resampledData(:,validFreqInds));
    power4fit = resampledData(:,validFreqInds);
    
    %Return the specgram to requested frequencies
    spec.freqs = spec.freqs(validFreqInds);
    spec.amp = spec.amp(:,validFreqInds);
    spec.data = spec.data(:,validFreqInds);
    
    spec.IRASAsmooth = power4fit;
    
    clear resampledData
else
   power4fit = spec.amp;
end

%% Fit the slope of the power spectrogram
rsq = zeros(size(spec.timestamps));
s = zeros(length(spec.timestamps),2);
yresid = zeros(length(spec.timestamps),length(spec.freqs));
for tt = 1:length(spec.timestamps)
    %Fit the line
    x = log10(spec.freqs);  y=power4fit(tt,:);
    s(tt,:) = polyfit(x,y,1);
    %Calculate the residuals (from the full PS)
    yfit =  s(tt,1) * x + s(tt,2);
    yresid(tt,:) = spec.amp(tt,:) - yfit; %residual between "raw" PS, not IRASA-smoothed
    %Calculate the rsquared value
    SSresid = sum(yresid(tt,:).^2);
    SStotal = (length(y)-1) * var(y);
    rsq(tt) = 1 - SSresid/SStotal;
end


%% Output Structure
specslope.data = s(:,1);
specslope.intercept = s(:,2);
specslope.timestamps = spec.timestamps;
specslope.specgram = spec.amp;
specslope.samplingRate = 1./dt;

specslope.detectionparms.winsize = winsize;
specslope.detectionparms.frange = frange;

specslope.rsq = rsq';
specslope.resid = yresid;
specslope.freqs = spec.freqs;

specslope.channels = lfp.channels;

if saveMat
    save(savename,'specslope','spec');
end

%% Figure
if SHOWFIG
    
    bigsamplewin = bz_RandomWindowInIntervals(spec.timestamps([1 end]),30);

   %hist(specslope.data,10)
   specmean.all = mean(spec.amp,1);
   slopebinIDs = discretize(specslope.data,linspace(min(specslope.data),max(specslope.data),6));
   for bb = 1:length(unique(slopebinIDs))
        specmean.bins(bb,:) = mean(spec.amp(slopebinIDs==bb,:),1);
        
        exwin(bb,:) = spec.timestamps(randsample(find(slopebinIDs==bb),1))+(winsize.*[-0.5 0.5]);
   end
figure
    subplot(4,1,1)
        imagesc(specslope.timestamps,log2(specslope.freqs),specslope.specgram')
        hold on
        plot(specslope.timestamps,bz_NormToRange(specslope.data),'w','linewidth',2)
        title('Spectrogram and PSS')
        LogScale('y',2)
        ylabel('f (Hz)')
        axis xy
        SpecColorRange( spec.amp );
        %caxis([0 1])
        xlim(bigsamplewin)
        bz_ScaleBar('s')
    subplot(4,1,2)
        imagesc(spec.timestamps,log2(spec.freqs),specslope.resid')
        hold on
        plot(lfp.timestamps,bz_NormToRange(single(lfp.data)),'w','linewidth',1)
        LogScale('y',2)
        ylabel('f (Hz)')
        title('Residuals and raw LFP')
        axis xy
        caxis([0 1])
        xlim(bigsamplewin)
        bz_ScaleBar('s')
%     subplot(8,1,3)
%         plot(lfp.timestamps,lfp.data,'k')
%         axis tight
%         box off
%         xlim(bigsamplewin)
%         lfprange = get(gca,'ylim');
%         set(gca,'XTickLabel',[])
%         ylabel('LFP')
         
    subplot(6,2,7)
        hist(specslope.data,10)
        box off
        xlabel('PSS')
        
        
    subplot(3,3,7)
        plot(log2(spec.freqs),specmean.all,'k','linewidth',2)
        hold on
        plot(log2(spec.freqs),specmean.bins,'k')
        LogScale('x',2)
        axis tight
        box off
        
    for bb = 1:2:5    
	subplot(6,2,12-(bb-1))
        plot(lfp.timestamps,lfp.data,'k')
        axis tight
        lfprange = get(gca,'ylim');
        box off
        xlim(exwin(bb,:)');ylim(lfprange)
        set(gca,'XTickLabel',[])
        set(gca,'YTick',[])
        ylabel('LFP')
    end
        bz_ScaleBar('s')

if saveMat
    figfolder = [basePath,filesep,'DetectionFigures'];
    NiceSave('PowerSpectrumSlope',figfolder,baseName)
end
        
end


end

