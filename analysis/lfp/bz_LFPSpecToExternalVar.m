function [ specvarcorr,specbyvar ] = bz_LFPSpecToExternalVar( LFP,extvar,specparms,figparms )
%[ specvarcorr,specbyvar] = bz_LFPSpecToExternalVar(LFP,extvar,specparms,figparms) 
%calculates the relationship of the LFP spectrum to an external 
%(e.g. behavioral) variable.
%Note: spectrogram are in units of dB (log(abs(spec)))
%
%INPUT
%   LFP         LFP structure in buzcode format:
%               a structure with fields lfp.data, lfp.timestamps
%   extvar      behavor structure in buzcode format:
%               a structure with fields extvar.data, extvar.timestamps
%               timestamps are in seconds, and aligned
%   specparms   structure of parameters for the spectrogram
%       .frange     [min max] frequency
%       .nfreqs     number of frequencies 
%       .space      spacing between freqs, 'lin' or 'log'
%       .specnorm   normalization for spectral power,
%                   options: 'mean','logmean','log' (default: log)
%       .numvarbins number of bins for your external variable
%       .varnorm    normalization for the variable,
%                   options: 'percentile', 'none' (default: 'none')
%       .type       'wavelet' or 'FFT'      (default: 'wavelet')
%      -if type: 'wavelet'-
%       .ncyc       number of wavelet cycles (recommended: ~5)
%      -if type: 'FFT'-
%       .winsize (s)
%       .overlap
%   figparms    (optional) parameters for the figure
%       .plotname   nametag for the saved figure
%       .figfolder  folder to save the figure in
%       .baseName   baseName of the recording for the figure
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
%
%DLevenstein 2017 (beta, feel free to improve...)
%TO DO IDEAS
%   -other spectrum types: FFT, MT
%   -default specparms
%% Calculate the spectrogram - FFT or WVLT
if ~isfield(specparms,'type')
    specparms.type = 'wavelet';
end

if ~isfield(specparms,'specnorm')
    specparms.specnorm = 'log';
end

switch specparms.type
    case 'wavelet'
        %Calcualte the Wavelet Transform
        [specvarcorr.freqs,~,spec] = WaveSpec(single(LFP.data),...
            specparms.frange,specparms.nfreqs,specparms.ncyc,...
            1/LFP.samplingRate,specparms.space);
        spectimestamps = LFP.timestamps; %Wavelet timestamp are same as LFP        

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
        spectimestamps = spectimestamps'+LFP.timestamps(1); %account for any time offset
        
        specbyvar.freqs = specvarcorr.freqs;
end

%Normalize the spectrogram
switch specparms.specnorm 
    case 'log'
        spec = log10(abs(spec));
    case 'mean'
        spec = (abs(spec));
        spec = bsxfun(@(X,Y) X./Y,spec,mean(spec,2));
    case 'logmean'
        spec = (abs(spec));
        spec = bsxfun(@(X,Y) X./Y,spec,mean(spec,2));
        spec = log10(spec);
end

%% Interpolate extvar to each LFP spectrogram timepoint
if ~isfield(specparms,'varnorm') %This is hacky
    specparms.varnorm ='none';
end

switch specparms.varnorm
    case 'percentile'
        vardata = NormToInt(extvar.data,'percentile');
    case 'none'
        vardata = extvar.data;
    otherwise
        warning('Your variable normalization input is wrong...')
end

%Get nearest value of extvar (behavior) for each point in the spectrogram
interpvar_LFP = interp1(extvar.timestamps,vardata,spectimestamps);

%% Binned Spectrogram
binedges = linspace(0,1,specparms.numvarbins+1);
specbyvar.varbins = binedges(1:specparms.numvarbins)+0.5.*diff(binedges([1 2]));
specbyvar.mean = zeros(specparms.nfreqs,specparms.numvarbins);
specbyvar.std = zeros(specparms.nfreqs,specparms.numvarbins);
for bb = 1:specparms.numvarbins
    inbinidx = interpvar_LFP>=binedges(bb) & interpvar_LFP<=binedges(bb+1);
    inbinspec = spec(:,inbinidx);
    specbyvar.mean(:,bb) = mean(inbinspec,2);
    specbyvar.std(:,bb) = std(inbinspec,[],2);
end

%% power-signal correlation for each frequency
sigthresh = 0.05; %Significance threshold (pvalue)
[specvarcorr.corr,specvarcorr.pval] = ...
    corr(interpvar_LFP(~isnan(interpvar_LFP)),...
    spec(:,~isnan(interpvar_LFP))','type','spearman');
corrsig = specvarcorr.pval<=sigthresh;

%% Figure
if exist('figparms','var')  %This whole figure thing can be better.
varcolormap=makeColorMap([0 0 0],[0.8 0 0],specparms.numvarbins);

figure
subplot(2,1,1)
    plot(extvar.timestamps,2*extvar.data-2,'k') %assumes 0-1 normalization... to fix later
    hold on
    imagesc(spectimestamps,log2(specvarcorr.freqs),spec)
    axis tight
    %hold on
    axis xy
    LogScale('y',2)
    switch specparms.specnorm
        case 'log'
            SpecColorRange(spec);
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
            ColorbarWithAxis([0.6 1.4],'Power (mean^1)');
        case 'logmean'
            ColorbarWithAxis([-0.1 0.1],'Power log10(mean^-^1)');
    end
    LogScale('y',2)
    ylabel('f (Hz)')
    
subplot(6,3,18)
    switch specparms.varnorm
        case 'percentile'
            plot(sort(vardata),sort(extvar.data),'r')
        case 'none'
            hist(interpvar_LFP(~isnan(interpvar_LFP)),specbyvar.varbins)
    end
    colorbar
    xlim(binedges([1 end]))
    xlabel('Behavior')
    set(gca,'ytick',[])
    
NiceSave(['SpecByVar',figparms.plotname],figparms.figfolder,figparms.baseName)
end

end

