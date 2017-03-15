function [phasedistros,phasebins,phasestats,h] = bz_PhaseModulation(varargin)
% USAGE
%[phasedistros,phasebins,phasestats,h] = bz_PhaseModulation(spikes,lfp,passband,intervals,samplingRate,method,plotting)
% 
% INPUTS
% spikes        - cell array where each element is a vector of
%                   spiketimes for each cell (time in seconds)
%
% lfp           - lfp values (1250kHz default), must be a single vector 
%
% passband      - frequency range for phase modulation [lowHz highHz] form
%
% intervals     - (optional) may specify timespans over which to calculate 
%               phase modulation.  Formats accepted: tstoolbox intervalSet
%               or a 2column matrix of [starts stops] in seconds
%
% samplingRate  - (optional) specifies lfp sampling frequency default:1250
%
% method        - (optional) method selection for how to generate phase, 
%               possibilties are: 'hilbert' (default) or 'wavelet'
%
% plotting      - (optional) 1 if want to plot, 0 if not. Default:1
%
%
% OUTPUTS
% phasedistros - Spike distribution perecentages for each cell in each bin
%                 specified by phasebins
% phasebins    - 180 bins spanning from 0 to 2pi
% phasestats   - ncellsx1 structure array with following (via
%                 CircularDistribution.m from FMAToolbox)
%                    phasestats.m        mean angle
%                    phasestats.mode     distribution mode
%                    phasestats.k        concentration
%                    phasestats.p        p-value for Rayleigh test
% h            - handles collection for all plotted figures
% 
%
%
% Calculates distribution of spikes over various phases from a specified
% cycle of an lfp vector.   Phase 0 means peak of lfp wave.
%
% Brendon Watson 2015
% edited by david tingley, 2017

%% defaults
p = inputParser;
addRequired(p,'spikes',@iscell);
addRequired(p,'lfp',@isnumeric);
addRequired(p,'passband',@isnumeric)

addParameter(p,'intervals',[0 inf],@isnumeric)
addParameter(p,'samplingRate',1250,@isnumeric)
addParameter(p,'method','hilbert',@isstr)
addParameter(p,'plotting',0,@isnumeric)
addParameter(p,'numBins',[180],@isnumeric)
% addParameter(p,'threshold',0,@isnumeric)

parse(p,varargin{:})

spikes = p.Results.spikes;
lfp = p.Results.lfp;
passband = p.Results.passband;

intervals = p.Results.intervals; % interval(s) over which to calculate
samplingRate = p.Results.samplingRate; % sampling rate of continuous signal (LFP)
method = p.Results.method; 
plotting = p.Results.plotting;
numBins = p.Results.numBins;



%% Get phase for every time point in LFP
switch lower(method)
    case ('hilbert')
        [b a] = butter(4,[passband(1)/(samplingRate/2) passband(2)/(samplingRate/2)],'bandpass');
%         [b a] = cheby2(4,20,passband/(samplingRate/2));
        filt = FiltFiltM(b,a,lfp);
        hilb = hilbert(filt);
        lfpphase = mod(angle(hilb),2*pi);
        clear fil
    case ('wavelet')% Use Wavelet transform to calulate the signal phases
        [wave,f,t,coh,wphases,raw,coi,scale,priod,scalef]=getWavelet(lfp(:,2),samplingRate,passband(1),passband(2),8,0);
        [~,mIdx]=max(wave);%get index max power for each timepiont
        pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
        lfpphase=wphases(pIdx);%get phase of max amplitude wave at each timepoint
        lfpphases = mod(lfpphases,2*pi);%covert to 0-2pi rather than -pi:pi
%     case ('peaks')
        % not yet coded
        % filter, smooth, diff = 0, diffdiff = negative
end
clear lfp

%% Get phases for each spike for each cell
h = [];
% cum_spkphases = [];
spkphases = cell(1,length(spikes));
for a = 1:length(spikes)
    
    bools = InIntervals(spikes{a},intervals);
    s =spikes{a}(bools);
%     s = spikes{a};
    if isempty(s)
        phasedistros(:,a) = zeros(numBins,1);
        phasestats.m(a) = nan;
        phasestats.r(a) = nan;
        phasestats.k(a) = nan;
        phasestats.p(a) = nan;
        phasestats.mode(a) = nan;
        spkphases{a} = nan;
    else
        spkphases{a} = lfpphase(round(s*samplingRate));

%         cum_spkphases = vertcat(cum_spkphases, spkphases{a});

    %% Gather binned counts and stats (incl Rayleigh Test)
        [phasedistros(:,a),phasebins,ps]=CircularDistribution(spkphases{a},'nBins',numBins);
        phasestats.m(a) = mod(ps.m,2*pi);
        phasestats.r(a) = ps.r;
        phasestats.k(a) = ps.k;
        phasestats.p(a) = ps.p;
        phasestats.mode(a) = ps.mode;

    %% plotting    
        if plotting
            h(end+1) = figure;
            hax = subplot(1,2,1); 
            rose(spkphases{a})
            title(hax,['Cell #' num2str(a) '. Rayleigh p = ' num2str(phasestats.p(a)) '.'])

            hax = subplot(1,2,2); 
            bar(phasebins*180/pi,phasedistros(:,a))
            xlim([0 360])
            set(hax,'XTick',[0 90 180 270 360]) 
            hold on;
            plot([0:360],cos(pi/180*[0:360])*0.05*max(phasedistros(:,a))+0.95*max(phasedistros(:,a)),'color',[.7 .7 .7])
            set(h(end),'name',['PhaseModPlotsForCell' num2str(a)]);
        end
    end
end

%% Cumulative effect across all spikes from all cells... not saving these stats for now
% phasebins=[];
% if length(cum_spkphases) > 10
%     [cpd,phasebins,cps]=CircularDistribution(cum_spkphases,'nBins',180);
%     cRp = cps.p;
% 
%     if plotting
%         h(end+1) = figure;
%         hax = subplot(1,2,1); 
%         rose(cum_spkphases)
%         title(hax,['All Spikes/Cells Accumulated. Rayleigh p = ' num2str(cps.p) '.'])
% 
%         hax = subplot(1,2,2); 
%         bar(phasebins*180/pi,cpd)
%         xlim([0 360])
%         set(hax,'XTick',[0 90 180 270 360]) 
%         hold on;
%         plot([0:360],cos(pi/180*[0:360])*0.05*max(cpd)+0.95*max(cpd),'color',[.7 .7 .7])
%         set(h(end),'name',['PhaseModPlotsForAllCells']);
%     end
% end

