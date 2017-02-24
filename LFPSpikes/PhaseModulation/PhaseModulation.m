function [phasedistros,phasebins,phasestats,h] = PhaseModulation(spikes,lfp,passband,intervals,lfpfreq,method,plotting)
% [phasedistros,phasebins,phasestats,h] = PhaseModulation(spikes,lfp,passband,intervals,lfpfreq,method,plotting)
% 
% Calculates distribution of spikes over various phases from a specified
% cycle of an lfp vector.   Phase 0 means peak of lfp wave.
%
% INPUTS
% spikes        - spike train, in seconds, a population or a single cell
%                   - may be a cell array where each element is a vector of
%                   spiketimes for each cell
%                   - may be tsdArray of multiple spike trains
%                   - may be a vector of spike times for a single cell
% lfp           - lfp values (1250kHz default), may be a single vector or
%               maybe in FMAtoolbox format of double vector with
%               times;values (only 2 columns not more)
% passband      - frequency range for phase modulation [lowHz highHz] form
%               (may be also other inputs valid for FilterLFP.m (FMAToolbox)
% intervals     - (optional) may specify timespans over which to calculate 
%               phase modulation.  Formats accepted: tstoolbox intervalSet
%               or a 2column matrix of [starts stops] in seconds
% lfpfreq       - (optional) specifies lfp sampling frequency default:1250
% method        - (optional) method selection for how to generate phase, 
%               possibilties are: 'hilbert' (default) or 'wavelet'
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
% Brendon Watson 2015



%% defaults
% if no intervals set to keep all recording
if ~exist('intervals','var')
    intervals = [0 Inf];
end
% get lfpfreq, if no input, assume 1250Hz
if ~exist('lfpfreq','var')
    lfpfreq = 1250;
end
% if no phase method specified, default is hilbert
if ~exist('method','var')
    method = 'hilbert';
end
% if no plotting preference specified, default is 1 (yes plot)
if ~exist('plotting','var')
    plotting = '1';
end

numbins = 180;


%% handle inputs, shape them as necessary
[spikes,lfp,intervals,lfpfreq] = handleInputsIn(spikes,lfp,intervals,lfpfreq);


%% Get phase for every time point in LFP
switch lower(method)
    case ('hilbert')
        fil=FilterLFP(lfp,'passband',passband); %from FMAToolbox
        hilb = hilbert(fil(:,2));
        lfpphase = mod(angle(hilb+pi),2*pi);
        clear fil
    case ('wavelet')% Use Wavelet transform to calulate the signal phases
        [wave,f,t,coh,wphases,raw,coi,scale,priod,scalef]=getWavelet(lfp(:,2),lfpfreq,passband(1),passband(2),8,0);
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
cum_spkphases = [];
spkphases = cell(1,length(spikes));
for a = 1:length(spikes)
    try
        bools = InIntervals(spikes{a},intervals);
    catch
        bools = InIntervalsBW(spikes{a},intervals);
    end
    s =spikes{a}(bools);
%     s = spikes{a};
    if isempty(s)
        phasedistros(:,a) = zeros(numbins,1);
        phasestats.m(a) = nan;
        phasestats.r(a) = nan;
        phasestats.k(a) = nan;
        phasestats.p(a) = nan;
        phasestats.mode(a) = nan;
        spkphases{a} = nan;
    else
        spkphases{a} = lfpphase(round(s));

        cum_spkphases = vertcat(cum_spkphases, spkphases{a});

    %% Gather binned counts and stats (incl Rayleigh Test)
        [phasedistros(:,a),phasebins,ps]=CircularDistribution(spkphases{a},'nBins',numbins);
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
[cpd,phasebins,cps]=CircularDistribution(cum_spkphases,'nBins',180);
cRp = cps.p;

if plotting
    h(end+1) = figure;
    hax = subplot(1,2,1); 
    rose(cum_spkphases)
    title(hax,['All Spikes/Cells Accumulated. Rayleigh p = ' num2str(cps.p) '.'])

    hax = subplot(1,2,2); 
    bar(phasebins*180/pi,cpd)
    xlim([0 360])
    set(hax,'XTick',[0 90 180 270 360]) 
    hold on;
    plot([0:360],cos(pi/180*[0:360])*0.05*max(cpd)+0.95*max(cpd),'color',[.7 .7 .7])
    set(h(end),'name',['PhaseModPlotsForAllCells']);
end


function [spikes,lfp,intervals,lfpfreq] = handleInputsIn(spikes,lfp,intervals,lfpfreq)
% Allow for multiple format options for certain inputs, convert them

% lfp AND lfpfreq setting
switch class(lfp)
    case {'double','single','int16'}
        lfp = double(lfp);
        if any(size(lfp) == 1);
            lfp = lfp(:);
            t = 1/lfpfreq:1/lfpfreq:length(lfp)/lfpfreq;
            t = t';
            lfp = cat(2,t,lfp);
            lfpfreq = 1/mode(diff(t));
        elseif any(size(lfp) == 2);
            %find which dim is 2
            sz = size(lfp);
            d = find(sz==2);
            if d == 1;
                t = lfp(1,:)';
                lfp = lfp(2,:)';
                lfp = cat(2,t,lfp);
            end                
            lfpfreq = 1/mode(diff(lfp(:,1)));
        end
        clear t
    case 'tsd'
        t = Range(lfp);
        t = t(:);
        lfp = Data(lfp);
        lfp = lfp(:);
        lfp = cat(2,t,lfp);
        lfpfreq = 1/mode(diff(t));
end

% convert spikes to a cell array of 1D spike timing vectors 
switch class(spikes)
    case 'tsdArray'
        S = spikes;
        spikes = {};
        for a = 1:length(S);
            t = Range(S{a},'s');
            spikes{a} = t(:)*lfpfreq;
        end
    case 'double'
        if numel(spikes) ~= length(spikes)
            error ('"spikes" input must be either 1D vector of spiketimes, tsdArray or cellarray of 1D vectors')
        else
            S = spikes;
            spikes = cell(1); % edited by Luke, 2015-08-10
            spikes{1} = S(:)*lfpfreq;
        end
    case 'cell'
        for a = 1:length(spikes)
            spikes{a} = spikes{a}(:)*lfpfreq;
        end
    otherwise
    error ('"spikes" input must be either 1D vector of spiketimes, tsdArray or cellarray of 1D vectors')
end


% get intervals if any and convert them to a pair of columns indicationg
switch class(intervals)
    case 'intervalSet'%just so reader knows this is an option
        I = intervals;
        intervals = [Start(I) End(I)];
    case 'double'
        if size(intervals,2) ~= 2
            error ('"intervals" input is optional but if entered it must be either an n x 2 vector of [start,stop] pairs or an intervalSet')
        end
    otherwise
            error ('"intervals" input is optional but if entered it must be either an n x 2 vector of [start,stop] pairs or an intervalSet')
end
intervals = intervals * lfpfreq;%convert all to timeframe of lfp


function [status,interval,index] = InIntervalsBW(times,ints)
% Slower version of FMAToolbox function InIntervals, but without dependence
% on compiled code which sometimes fails.
% 
% InIntervals - Test which values fall in a list of intervals.
%
%  USAGE
%
%    [status,interval,index] = InIntervals(values,intervals,<options>)
%
%    values         values to test (these need not be ordered)
%    intervals      list of (start,stop) pairs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'verbose'     display information about ongoing processing
%                   (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    status         logical indices (1 = belongs to one of the intervals,
%                   0 = belongs to none)
%    interval       for each value, the index of the interval to which
%                   it belongs (0 = none)
%    index          for each value, its index in the interval to which
%                   it belongs (0 = none)
%
% Brendon Watson 2015

allintids = [];
allidxs = [];

times = times(:);
for a = 1:size(ints,1)
    tbool = times>=ints(a,1) & times<=ints(a,2);
    tintids = a*tbool;
    tidxs = cumsum(tbool);
    tidxs(~tbool) = 0;

    allintids = cat(2,allintids,tintids);
    allidxs = cat(2,allidxs,tidxs);
end

status = logical(sum(allintids,2));
interval = sum(allintids,2);
index = sum(allidxs,2);
