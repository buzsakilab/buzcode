function [phaseamplitude] = bz_PhaseAmplitudeDist(lfp,phaserange,amprange,varargin)
% [phaseamplitudemap,ampfreqs,phasecenters] = bz_PhaseAmplitudeDist(lfp,phaserange,amprange)
%This function calculates the mean amplitude of multiple higher frequency bands at
%a the phase for a single lower frequency band signal
%
%%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%    phaserange     [min max] frequency to filter for the single phase signal
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
%   phaseamplitude   phase by frequency matrix. value is mean amplitude 
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
%Last Updated: 31/10/2019 (EFO)
%DLevenstein

%
%
%DLevenstein 2016
%Updated by Eliezyer de Oliveira (EFO)/2019

%% parsing inputs

if ~bz_isLFP(lfp)
    error('Not following the lfp structure required, see documentation')
end


p = inputParser;
addParameter(p,'filterType','fir1',@ischar);
addParameter(p,'filterOrder',4,@isnumeric);
addParameter(p,'numBins',50,@isnumeric);
addParameter(p,'intervals',[0 inf],@isnumeric);
addParameter(p,'nfreqs',100,@isnumeric);
addParameter(p,'ncyc',7,@isnumeric);
addParameter(p,'makePlot',true,@islogical);

parse(p,varargin{:});
makePlot = p.Results.makePlot;
filterType = p.Results.filterType;
filterOrder = p.Results.filterOrder;
numBins = p.Results.numBins;
nfreqs = p.Results.nfreqs;
ncyc = p.Results.ncyc;
intervals = p.Results.intervals;

%% Deal with input types
% sf_LFP = lfp.samplingRate;


%% Filter LFP for the phase
filtered_phase = bz_Filter(lfp,'passband',phaserange,'filter',filterType,'order',filterOrder,'intervals',intervals);

%% Get LFP, Phase in intervals
%edgebuffer = 1; %s
%edgebuffer_si = edgebuffer.*sf_LFP;
%edgebuffer = edgebuffer.*[1 1];
% LFP_int = IsolateEpochs2(LFP.data,int,edgebuffer,sf_LFP);
% LFPphase_int = IsolateEpochs2(filtered.phase,int,edgebuffer,sf_LFP);

%% Wavelet Transform LFP in intervals
wavespec_amp = bz_WaveSpec(lfp,'frange',amprange,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',intervals);
%[ampfreqs,~,spec_int]
%spec_int = cellfun(@(X) abs(X),spec_int,'UniformOutput',false);
wavespec_amp.data = log10(abs(wavespec_amp.data)); %log scaled for normality
wavespec_amp.mean = mean(wavespec_amp.data,1);
%Remove Buffers (this should be done via 'intervals' in bz_WaveSpec
%spec_int = cellfun(@(X) X(:,edgebuffer_si:end-edgebuffer_si),spec_int,'UniformOutput',false);
%LFPphase_int = cellfun(@(X) X(edgebuffer_si:end-edgebuffer_si),LFPphase_int,'UniformOutput',false);

%% Bin phase and power
phasebins = linspace(-pi,pi,numBins+1);
phasecenters = phasebins(1:end-1)+(diff(phasebins)/2);

[~,~,phaseall] = histcounts(filtered_phase.phase,phasebins);

phaseamplitudemap = zeros(numBins,nfreqs);
for bb = 1:numBins
    phaseamplitudemap(bb,:) = mean(wavespec_amp.data(phaseall==bb,:),1)./wavespec_amp.mean;
end
 
ampfreqs = wavespec_amp.freqs;

phaseamplitude.map = phaseamplitudemap;
phaseamplitude.phasecenters = phasecenters;
phaseamplitude.amp_freq = ampfreqs;
phaseamplitude.phase_range = num2str(mean(phaserange));
phaseamplitude.params.filterType  = filterType;
phaseamplitude.params.filterOrder = filterOrder;

%% Plot
if makePlot
    figure
        imagesc(phaseamplitude.phasecenters,log2(phaseamplitude.amp_freq),((phaseamplitude.map))')
        hold on
        imagesc(phaseamplitude.phasecenters+2*pi,log2(phaseamplitude.amp_freq),((phaseamplitude.map))')
        plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(phaseamplitude.amp_freq(end/2)),'k')
        xlabel(['Phase (',num2str(phaseamplitude.phase_range),' Hz)']);
        ylabel('Amplitude Frequency (Hz)')
        LogScale('y',2)
        xlim([phaseamplitude.phasecenters(1) phaseamplitude.phasecenters(end)+2*pi]);
        colorbar
        axis xy
end

end

