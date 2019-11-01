function [comod] = bz_CrossFreqMod(lfp,phaserange,amprange,varargin)
% [comod] = bz_CrossFreqMod(lfp,phaserange,amprange,flagPlot)
%
%
%This function calculates the modulation index of phase-amplitude between
%phaserange (lower frequencies) to amplitude range (higher frequencies).
%It can really take a long time if you do very small steps of frequency,
%due to wavelet processing each frequency at a time.
%
%%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    phaserange     [min:steps:max] array of frequencies to filter for the
%                   phase signal
%    amprange       [min:stepsmax] array of frequencies range for wavelets
%                   for the power signal
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     phaseCh      channel to compute phase. If empty take first channel
%     ampChans     channels to compute amplitude. If empty take first channel
%     makePlot      default true
%    =========================================================================
%
%OUTPUT
%   comod               phase-frequency versus amplitude-frequency
%                       comodulogram matrix in modulation index units
%
%Dependencies
%   bz_Filter
%   bz_WaveSpec
%
%   AntonioFR, 2/2019

%% Parse inputs

p = inputParser;
addParameter(p,'phaseCh',lfp.channels(1),@isnumeric);
addParameter(p,'ampChans',lfp.channels(1),@isnumeric);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'makePlot',true,@islogical);

parse(p,varargin{:});
phaseCh = p.Results.phaseCh;
ampChans = p.Results.ampChans;
samplingRate = p.Results.samplingRate;
makePlot = p.Results.makePlot;

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

%% Filter LFP for phase
nfreqs = length(amprange);

for bnd = 1:length(phaserange)-1
    filtered_phase(bnd,:) = bz_Filter(lfp,'passband',phaserange(bnd:bnd+1),'filter','fir1','channels',phaseCh);
end

%% Wavelet Transform LFP in intervals
for ch = 1:length(ampChans)
    comod = zeros(length(amprange)-1,length(filtered_phase),length(ampChans));
    
    lfpCh = lfp; % this should be remove when bz_WaveSpec can take 'channels' input
    lfpCh.data = lfp.data(:,find(lfp.channels == ampChans(ch)));
    
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
            comod(apr,idx,ch) = sum(phaseAmp.*log(phaseAmp./(ones(numbins,size(phaseAmp,2))/numbins)))/log(numbins);
        end

    end
    ampfreqs = wavespec_amp.freqs;
    clear lfpCh
end


%% Plot
if makePlot
    
    figure
    for ch = 1:length(ampChans)
        subplot(1,length(ampChans),ch);
        contourf(phaserange(2:end),amprange(2:end),abs(comod(:,:,ch)),20,'LineColor','none');
        colorbar ('SouthOutside'); colormap jet;
        %title(signals)
        if ch > 1
            set(gca,'YTick',[]);
        end
    end
    
end

