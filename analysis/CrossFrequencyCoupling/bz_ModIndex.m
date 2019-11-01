function [comod] = bz_ModIndex(lfp,phaserange,amprange,flagPlot)
% [comod] = bz_ModIndex(lfp,phaserange,amprange,flagPlot)
%This function calculates the modulation index of phase-amplitude between
%phaserange (lower frequencies) to amplitude range (higher frequencies).
%It can really take a long time if you do very small steps of frequency,
%due to wavelet processing each frequency at a time.
%
%%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%    phaserange     [min:steps:max] array of frequencies to filter for the
%                   phase signal
%    amprange       [min:stepsmax] array of frequencies range for wavelets
%                   for the power signal
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
%   comod               Modulation index matrix between phase frequency and
%                       amplitude
%
%Dependencies
%   bz_Filter
%   bz_WaveSpec
%
% Implemented by Eliezyer de Oliveira, 2018
% Last Update: 08/06/2018
%
% I couldn't found a better/faster way to do this processing, I could run
% all the frequencies at once but if the recording is too long it might not
% have enough memory to hold all the data.




nfreqs = length(amprange);

%% Filter LFP for phase
parfor bnd = 1:length(phaserange)-1
    filtered_phase(bnd,:) = bz_Filter(lfp,'passband',phaserange(bnd:bnd+1),'filter','fir1');
end

%% Wavelet Transform LFP in intervals
comod = zeros(length(amprange)-1,length(filtered_phase));
for apr = 1:length(amprange)-1
    wavespec_amp = bz_WaveSpec(lfp,'frange',[amprange(apr) amprange(apr+1)],'nfreqs',1);
    
    wavespec_amp.data = abs(wavespec_amp.data);
    %% Bin phase and power
    numbins = 50;
    phasebins = linspace(-pi,pi,numbins+1);
    phasecenters = phasebins(1:end-1)+(phasebins(2)-phasebins(1));
    
    for idx = 1:length(filtered_phase)
        [phasedist,~,phaseall] = histcounts(filtered_phase(idx).phase,phasebins);
        
        phaseAmp = zeros(numbins,1);
        for bb = 1:numbins
            phaseAmp(bb) = mean(wavespec_amp.data(phaseall==bb),1);
        end
        
        phaseAmp = phaseAmp./sum(phaseAmp,1);
        comod(apr,idx) = sum(phaseAmp.*log(phaseAmp./(ones(numbins,size(phaseAmp,2))/numbins)))/log(numbins);
    end
    
end

ampfreqs = wavespec_amp.freqs;
%% Plot
if flagPlot
    figure
    imagesc(phaserange,log2(ampfreqs),comod);
    colormap jet
    hold on
    xlabel('Frequency phase');
    ylabel('Frequency amplitude')
    LogScale('y',2)
    colorbar
    axis xy
    
end

