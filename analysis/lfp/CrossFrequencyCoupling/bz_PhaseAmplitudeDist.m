function [phaseamplitudemap,ampfreqs,phasecenters] = bz_PhaseAmplitudeDist(lfp,phaserange,amprange)
% [phaseamplitudemap,ampfreqs,phasecenters] = bz_PhaseAmplitudeDist(lfp,phaserange,amprange)
%This function calculates the mean amplitude of higher frequency bands at
%a the phase for a given lower frequency band signal
%
%%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%    phaserange     [min max] frequency to filter for the phase signal
%    amprange       [min max] frequency range for wavelets for the power
%                   signal
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       not implemented yet.... 'nfreqs' for amplitude signal
%                               'ncyc' for wavelet parameters
%                               'interval' interval on which to calculate
%                               other filter parameters
%                               'phasesignal' option for other signal to be
%                                               used as phase information
%    =========================================================================
%
%OUTPUT
%   phaseamplitudemap   phase by frequency matrix. value is mean amplitude 
%                       of that frequency at that phase 
%                       (relative to the mean over the entire signal)
%   ampfreqs            amplitude frequencies
%   phasecenters        phase bins
%
%
%Dependencies
%   bz_Filter
%   bz_WaveSpec
%
%Last Updated: 10/9/15
%DLevenstein

%
%
%DLevenstein 2016
%TO DO: intervals support
%% DEV (these should be varagins)
nfreqs = 100;
ncyc = 7;

%% Deal with input types
sf_LFP = lfp.samplingRate;


%% Filter LFP for the phase
filtered_phase = bz_Filter(lfp,'passband',phaserange,'filter','fir1');

%% Get LFP, Phase in intervals
%edgebuffer = 1; %s
%edgebuffer_si = edgebuffer.*sf_LFP;
%edgebuffer = edgebuffer.*[1 1];
% LFP_int = IsolateEpochs2(LFP.data,int,edgebuffer,sf_LFP);
% LFPphase_int = IsolateEpochs2(filtered.phase,int,edgebuffer,sf_LFP);

%% Wavelet Transform LFP in intervals
wavespec_amp = bz_WaveSpec(lfp,'frange',amprange);
%[ampfreqs,~,spec_int]
%spec_int = cellfun(@(X) abs(X),spec_int,'UniformOutput',false);
wavespec_amp.data = log10(abs(wavespec_amp.data)); %log scaled for normality
wavespec_amp.mean = mean(wavespec_amp.data,1);
%Remove Buffers (this should be done via 'intervals' in bz_WaveSpec
%spec_int = cellfun(@(X) X(:,edgebuffer_si:end-edgebuffer_si),spec_int,'UniformOutput',false);
%LFPphase_int = cellfun(@(X) X(edgebuffer_si:end-edgebuffer_si),LFPphase_int,'UniformOutput',false);

%% Bin phase and power
numbins = 50;
phasebins = linspace(-pi,pi,numbins+1);
phasecenters = phasebins(1:end-1)+(phasebins(2)-phasebins(1));

[phasedist,~,phaseall] = histcounts(filtered_phase.phase,phasebins);


phaseamplitudemap = zeros(numbins,nfreqs);
for bb = 1:numbins
    phaseamplitudemap(bb,:) = mean(wavespec_amp.data(phaseall==bb,:),1)./wavespec_amp.mean;
end
 
ampfreqs = wavespec_amp.freqs;
%% Plot
figure
    imagesc(phasecenters,log2(ampfreqs),((phaseamplitudemap))')
    hold on
    imagesc(phasecenters+2*pi,log2(ampfreqs),((phaseamplitudemap))')
    plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(ampfreqs(end/2)),'k')
    xlabel(['Phase (',num2str(phaserange),'Hz)']);
    ylabel('Frequency (Hz')
    LogScale('y',2)
    xlim([phasecenters(1) phasecenters(end)+2*pi]);
    colorbar
    axis xy

end

