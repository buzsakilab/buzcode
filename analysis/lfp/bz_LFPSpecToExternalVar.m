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
%       .ncyc       number of wavelet cycles (recommended: ~5)
%       .space      spacing between freqs, 'lin' or 'log'
%       .numvarbins number of bins for your external variable
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
%% Calcualte the Wavelet Transform
[specvarcorr.freqs,t,spec] = WaveSpec(single(LFP.data),...
    specparms.frange,specparms.nfreqs,specparms.ncyc,...
    1/LFP.samplingRate,specparms.space);
spec = log10(abs(spec));

specbyvar.freqs = specvarcorr.freqs;
%% Interpolate time
%Get nearest value of extvar (behavior) for each point in the LFP
interpvar_LFP = interp1(extvar.timestamps,extvar.data,LFP.timestamps);

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

%% correlation for each frequency
sigthresh = 0.05;
[specvarcorr.corr,specvarcorr.pval] = corr(interpvar_LFP(~isnan(interpvar_LFP)),spec(:,~isnan(interpvar_LFP))','type','spearman');
corrsig = specvarcorr.pval<=sigthresh;
%% Figure
if exist('figparms','var')
varcolormap=makeColorMap([0 0 0],[0.8 0 0],specparms.numvarbins);

figure
subplot(2,1,1)
    plot(extvar.timestamps,extvar.data-1,'k')
    hold on
    imagesc(t,log2(specvarcorr.freqs),spec)
    axis tight
    %hold on
    axis xy
    LogScale('y',2)
    SpecColorRange(spec);
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
    %colorbar
    LogScale('y',2)
    ylabel('f (Hz)')
    
subplot(6,3,18)
    hist(interpvar_LFP(~isnan(interpvar_LFP)),specbyvar.varbins)
    xlim(binedges([1 end]))
    xlabel('Behavior')
    set(gca,'ytick',[])
    
NiceSave(['SpecByVar',figparms.plotname],figparms.figfolder,figparms.baseName)
end

end

