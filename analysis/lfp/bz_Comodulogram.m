function [ comod ] = bz_Comodulogram(lfp,specparms,figparms)
%[ comod ] = bz_Comodulogram(lfp,specparms,figparms) calculates the
%comodulogram (i.e. power-power correlation) for an lfp file.
%
%
%INPUT
%   LFP         LFP structure in buzcode format:
%               a structure with fields lfp.data, lfp.timestamps
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
%DLevenstein 2017
%TO DO 
% -update to work with multiple LFP channels...
% -use input parser to include default specparms
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
        [comod.freqs,~,spec] = bz_WaveSpec(single(lfp.data),...
            specparms.frange,specparms.nfreqs,specparms.ncyc,...
            1/lfp.samplingRate,specparms.space);
        spectimestamps = LFP.timestamps; %Wavelet timestamp are same as LFP        
        
	case 'FFT'
        %Calculate the frequences to use
        switch specparms.space
            case 'log'
                comod.freqs = logspace(log10(specparms.frange(1)),...
                    log10(specparms.frange(2)),specparms.nfreqs);
            case 'lin'
                comod.freqs = linspace(specparms.frange(1),...
                    specparms.frange(2),specparms.nfreqs);  
        end
        
        %Calculate the FFT spectrogram parameters - covert from s to sf
        winsize = specparms.winsize*lfp.samplingRate;
        noverlap = specparms.noverlap*lfp.samplingRate;
        %Calculate the FFT spectrogram
        [spec,~,spectimestamps] = spectrogram(single(lfp.data),...
            winsize,noverlap,comod.freqs,lfp.samplingRate);
        spectimestamps = spectimestamps'+lfp.timestamps(1); %account for any time offset

end

spec = log10(abs(spec)); %Log-transform power
%% Calculate the power-power correlations
comod.corrs = corr(spec','type','spearman');

%%
if exist('figparms','var')  %This whole figure thing can be better.
corrcolor= [makeColorMap([1 1 1],[0 0 0.8],[0 0 0]);...
    makeColorMap([0 0 0],[0.8 0 0],[1 1 1])];
figure
colormap(corrcolor)
    imagesc(log2(comod.freqs),log2(comod.freqs),comod.corrs)
    colorbar
    ColorbarWithAxis([-0.4 0.4],'Power-Power Correlation (rho)')
    LogScale('xy',2)
    xlabel('f (Hz)');ylabel('f (Hz)')
    
NiceSave(['Comodulogram',figparms.plotname],figparms.figfolder,figparms.baseName)
end

