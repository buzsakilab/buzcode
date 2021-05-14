function [ specvarcorr,specbyvar ] = bz_LFPSpecToExternalVar( LFP,extvar,varargin )
%[ specvarcorr,specbyvar] = bz_LFPSpecToExternalVar(LFP,extvar,<options>) 
%calculates the relationship of the LFP spectrum to an external 
%(e.g. behavioral) variable.
%Note: spectrogram are in units of dB (log(abs(spec)))
%
%INPUT
%   LFP         LFP structure in buzcode format:
%               a structure with fields lfp.data, lfp.timestamps
%   extvar      behavor structure in buzcode format:
%               a structure with fields extvar.data, extvar.timestamps
%               timestamps are in seconds, and aligned.
%               if behavior is a vector, it must be time-aligned with LFP
%
%   <options>
%   'specparms'   structure of parameters for the spectrogram. 
%                 use 'specparms','input' if spectrogram is provided instead
%                 of lfp
%       .frange     [min max] frequency  (default: [1 128])
%       .nfreqs     number of frequencies (default: 100)
%       .space      spacing between freqs, 'lin' or 'log' (default: log)
%       .specnorm   normalization for spectral power,
%                   options: 'mean','logmean','log' (default: log)
%       .type       'wavelet' or 'FFT'      (default: 'wavelet')
%      -if type: 'wavelet'-
%       .ncyc       number of wavelet cycles (recommended: ~5)
%      -if type: 'FFT'-
%       .winsize (s)
%       .noverlap
%   'figparms'      (optional) parameters for the figure
%       .plotname   nametag for the saved figure
%       .figfolder  folder to save the figure in
%       .baseName   baseName of the recording for the figure
%   'numvarbins'      number of bins for your external variable
%   'varlim'          [lower upper] bounds for variable bins (default: min/max)
%   'varnorm'         normalization for the variable,
%                   options: 'percentile', 'none' (default: 'none')
%   'ints'          intervals to restrict analysis to
%   'minX'          minimum number of timepoints needed to calculate mean (default: 25)
%
%OUTPUT
%   specvarcorr     correlation between power in each frequency and the var
%       .corr       correlation value (spearman) for each frequency
%       .pval       p value for the correlation for each frequency
%       .freqs      frequencies (Hz)
%   specbyvar       wavelet spectrogram w.r.t binned value of behavior
%       .mean       [nfreqs x nbins] mean spectrogram in each behavior bin 
%       .std        [nfreqs x nbins] std spectrogram in each behavior bin
%       .freqs      frequencies (Hz)
%       .varbins    behvior bin centers
%       .vardist    portion of time in each variable bin
%
%DLevenstein 2017 (beta, feel free to improve...)
%TO DO IDEAS
%   -other spectrum types: MT
%   -multi-channel lfp... just loop each channel/column, select 'channels'
%   -default specparms
%%
defaultspecparms.frange = [1 128];
defaultspecparms.nfreqs = 100;
defaultspecparms.space = 'log';
defaultspecparms.specnorm = 'log';
defaultspecparms.type = 'wavelet';
defaultspecparms.ncyc = 5;

p = inputParser;
addParameter(p,'specparms',defaultspecparms);
addParameter(p,'figparms',[]);
addParameter(p,'numvarbins',10);
addParameter(p,'varlim',[]);
addParameter(p,'varnorm','none');
addParameter(p,'ints',[0 Inf]);
addParameter(p,'minX',25);

parse(p,varargin{:})
specparms = p.Results.specparms;
figparms = p.Results.figparms;
numvarbins = p.Results.numvarbins;
varnorm = p.Results.varnorm;
varlim = p.Results.varlim;
ints = p.Results.ints;
minX = p.Results.minX;


%% Calculate the spectrogram - FFT or WVLT
if strcmp(specparms,'input')
    clear specparms
    if isstruct(LFP)
        spec = LFP.data;
        spectimestamps = LFP.timestamps;
        specparms.freqs = LFP.freqs;
    else
       spec = LFP;
       spectimestamps = [1:size(LFP,1)]';
       specparms.freqs = 1:size(LFP,2);
    end
    
    
    specvarcorr.freqs = specparms.freqs;
    specparms.nfreqs = length(specparms.freqs);
    specparms.specnorm = 'log';
else
    if ~isfield(specparms,'type')
        specparms.type = 'wavelet';
    end

    if ~isfield(specparms,'specnorm')
        specparms.specnorm = 'log';
    end

    switch specparms.type
        case 'wavelet'
            %Calcualte the Wavelet Transform
            [wavespec] = bz_WaveSpec(LFP,...
                'frange',specparms.frange,'nfreqs',specparms.nfreqs,...
                'ncyc',specparms.ncyc,'space',specparms.space);

            spectimestamps = wavespec.timestamps; %Wavelet timestamp are same as LFP        
            spec = wavespec.data;
            specvarcorr.freqs = wavespec.freqs;
            specbyvar.freqs = specvarcorr.freqs;

        case 'FFT'
            %Calculate the frequences to use
            switch specparms.space
                case 'log'
                    specvarcorr.freqs = logspace(log10(specparms.frange(1)),...
                        log10(specparms.frange(2)),specparms.nfreqs);
                case 'lin'
                    specvarcorr.freqs = linspace(specparms.frange(1),...
                        specparms.frange(2),specparms.nfreqs);  
            end

            %Calculate the FFT spectrogram parameters - covert from s to sf
            winsize = specparms.winsize*LFP.samplingRate;
            noverlap = specparms.noverlap*LFP.samplingRate; %Might want to calaulte overlap based on pupil...?
            %Calculate the FFT spectrogram
            [spec,~,spectimestamps] = spectrogram(single(LFP.data),...
                winsize,noverlap,specvarcorr.freqs,LFP.samplingRate);
            spec = spec';
            spectimestamps = spectimestamps'+LFP.timestamps(1); %account for any time offset

            specbyvar.freqs = specvarcorr.freqs;
    end

    %Normalize the spectrogram
    switch specparms.specnorm 
        case 'log'
            spec = log10(abs(spec));
        case 'mean'
            spec = (abs(spec));
            spec = bsxfun(@(X,Y) X./Y,spec,mean(spec,1));
        case 'logmean'
            spec = (abs(spec));
            spec = bsxfun(@(X,Y) X./Y,spec,mean(spec,1));
            spec = log10(spec);
    end

end

%% Variable - is it a vector
if ~isstruct(extvar)
    temp.data = extvar;
    temp.timestamps = spectimestamps;
    extvar = temp;
end
%% Interpolate extvar to each LFP spectrogram timepoint
switch varnorm
    case 'percentile'
        vardata = NormToInt(extvar.data,'percentile');
    case 'none'
        vardata = extvar.data;
    otherwise
        warning('Your variable normalization input is wrong...')
end

%Make spectrogram and extvar on the same time points
switch length(extvar.timestamps) >= length(spectimestamps)
    case true  %If fewer timepointes in the spectrogram, 
               %time to spec timepoints that are close to vardatatimepoints
        %Find closest variable time point for each spectrogram time point
        specatvartimes = interp1(extvar.timestamps,extvar.timestamps,spectimestamps,'nearest');
        specatvartimes = unique(specatvartimes(~isnan(specatvartimes)));
        %Remove any doubles/nans (i.e. intervals where there are no extvar timepoints
        vardata = interp1(extvar.timestamps,vardata,specatvartimes,'nearest');
        interpspec = interp1(spectimestamps,spec,specatvartimes,'nearest');
    case false %If fewer timepointes in extvar, interp time to extvar
        interpspec = interp1(spectimestamps,spec,extvar.timestamps,'nearest');
end
%Get nearest value of extvar (behavior) for each point in the spectrogram
%vardata = interp1(extvar.timestamps,vardata,spectimestamps);
%Get the nearest value of the spectrogram for each timepoint in extvar (behavior)
%interpspec = interp1(spectimestamps,spec',extvar.timestamps,'nearest');
%% REstrict to intervals

inints = InIntervals(specatvartimes,ints);
vardata = vardata(inints);
interpspec = interpspec(inints,:);
%% Binned Spectrogram
if isempty(varlim)
    varlim(1) = min(vardata); varlim(2)=max(vardata);
end
    
binedges = linspace(varlim(1),varlim(2),numvarbins+1);
specbyvar.varbins = binedges(1:numvarbins)+0.5.*diff(binedges([1 2]));
specbyvar.mean = zeros(specparms.nfreqs,numvarbins);
specbyvar.std = zeros(specparms.nfreqs,numvarbins);
for bb = 1:numvarbins
    inbinidx = vardata>=binedges(bb) & vardata<=binedges(bb+1);
    if sum(inbinidx)<minX
        specbyvar.mean(:,bb) = nan;
        specbyvar.std(:,bb) = nan;
        continue
    end
    inbinspec = interpspec(inbinidx,:);
    specbyvar.mean(:,bb) = mean(inbinspec,1)';
    specbyvar.std(:,bb) = std(inbinspec,[],1)';
end

specbyvar.vardist = hist(vardata(~isnan(vardata)),specbyvar.varbins);
specbyvar.vardist = specbyvar.vardist./sum(specbyvar.vardist);
%% power-signal correlation for each frequency
sigthresh = 0.05; %Significance threshold (pvalue)
% [specvarcorr.corr,specvarcorr.pval] = ...
%     corr(interpvar_LFP(~isnan(interpvar_LFP)),...
%     spec(:,~isnan(interpvar_LFP))','type','spearman');
[specvarcorr.corr,specvarcorr.pval] = ...
    corr(vardata,interpspec,'rows','complete','type','spearman');
corrsig = specvarcorr.pval<=sigthresh;

%% Figure
if ~isempty(figparms)  %This whole figure thing can be better.
varcolormap=makeColorMap([0 0 0],[0.8 0 0],numvarbins);

var4plot = (extvar.data-varlim(1))./(varlim(2)-varlim(1));

figure
subplot(2,1,1)
    plot(extvar.timestamps,2*var4plot-2,'k.') %assumes 0-1 normalization... to fix later
    hold on
    imagesc(spectimestamps,log2(specvarcorr.freqs),spec')
    axis tight
    %hold on
    axis xy
    LogScale('y',2)
    switch specparms.specnorm
        case 'log'
            try
                SpecColorRange(spec);
            end
        case 'mean'
            ColorbarWithAxis([0 4],'Power (mean^-^1)');
        case 'logmean'
            ColorbarWithAxis([-1 1],'Power log10(mean^-^1)');
    end
    xlabel('t (s)');ylabel('f (Hz)')
subplot(2,3,4)
    plot(log2(specvarcorr.freqs),specvarcorr.corr)
    hold on
    plot(log2(specvarcorr.freqs(corrsig)),specvarcorr.corr(corrsig),'.r')
    plot(log2(specvarcorr.freqs([1 end])),[0 0],'r--')
    LogScale('x',2)
    xlabel('Frequency (f)');ylabel('Pupil-Power Correlation (rho)')
    axis tight
    
subplot(2,3,5)
    set(gca,'colororder',varcolormap)
    hold all
    plot(log2(specvarcorr.freqs),specbyvar.mean)
    axis tight
    LogScale('x',2)
    xlabel('f (Hz)');ylabel('Power (dB)')
    
subplot(6,3,[12 15])
    imagesc(specbyvar.varbins,log2(specvarcorr.freqs),(specbyvar.mean))
    axis xy
    
    switch specparms.specnorm
        case 'log'
            colorbar
        case 'mean'
            ColorbarWithAxis([0.4 1.6],'Power (mean^1)');
            colormap(gca,'jet')
        case 'logmean'
            ColorbarWithAxis([-0.1 0.1],'Power log10(mean^-^1)');
    end
    LogScale('y',2)
    ylabel('f (Hz)')
    
subplot(6,3,18)
    switch varnorm
        case 'percentile'
            plot(sort(vardata),sort(extvar.data),'r')
        case 'none'
            hist(vardata(~isnan(vardata)),specbyvar.varbins)
    end
    colorbar
    xlim(binedges([1 end]))
    xlabel('Behavior')
    set(gca,'ytick',[])
    
    if isfield(figparms,'figfolder')
        NiceSave(['SpecByVar',figparms.plotname],figparms.figfolder,figparms.baseName)
    end
end
colormap jet
end



