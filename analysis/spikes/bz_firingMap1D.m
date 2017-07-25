function [firingMaps] = bz_firingMap1D(varargin)
% USAGE
% [rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes,behavior,lfp,tau)
%
% INPUTS
%
%   spikes  - buzcode format .cellinfo. struct with the following fields
%           .times 
%
%   behavior
%
%   lfp
%
%   tau
%
%
% OUTPUT
%  phaseMap - 7xM matrix columns are linearized position, trial
%  number, x position, y position, instananeous firing rate, theta phase.  M is the number of spikes
%  across all trials of the same type
%
%
%
%
%
%
% written by david tingley, 2017


% parse inputs
p=inputParser;
addRequired(p,'spikes',@isstruct);
addRequired(p,'behavior',@isstruct);
addRequired(p,'lfp',@isstruct);
addRequired(p,'tau',@isnumeric);
addParameter(p,'saveMat',false,@islogical);
parse(p,varargin{:});

spikes = p.Results.spikes;
behavior = p.Results.behavior;
lfp = p.Results.lfp;
tau = p.Results.tau;
saveMat = p.Results.saveMat;

% start processing
[b a] = butter(3,[6/625 9/625],'bandpass');
phases = (angle(hilbert(filtfilt(b,a,double(lfp.data(:,1))))));

for tt =1:length(unique(behavior.events.trialConditions))
    trials = find(behavior.events.trialConditions==tt);
    rateMap{tt} = zeros(length(spikes.times),length(trials),length(behavior.events.map{tt}.x));
    countMap{tt} = zeros(length(spikes.times),length(trials),length(behavior.events.map{tt}.x));
    occuMap{tt} = zeros(length(trials),length(behavior.events.map{tt}.x));
    if ~isempty(behavior.events.map{tt})
        for i =1:length(spikes.times)
            phaseMap{tt}{i} = [];
        end

        for i =1:length(spikes.times)
            for t=1:length(trials)
                if ~isempty(spikes.times{i})
                start = behavior.events.trialIntervals(trials(t),1);
                stop =  behavior.events.trialIntervals(trials(t),2);
                f = find(spikes.times{i}>start);
                ff = find(spikes.times{i}<stop);

                sp = intersect(f,ff);  % all spikes.times from a single trials
                sp_w_ends = [start spikes.times{i}(sp)' stop ];
                ISI = diff(sp_w_ends);
                for s=1:length(ISI)-1
                   inst_rate(s) = 1./mean(ISI(s:s+1)); 
                end 
                for s =1:length(sp)
                    [a b] = min(abs(spikes.times{i}(sp(s))-behavior.events.trials{trials(t)}.timestamps(:,1)));

                    x = behavior.events.trials{trials(t)}.x(b);
                    y = behavior.events.trials{trials(t)}.y(b);
                    [aa bb] = min(abs(behavior.events.map{tt}.x-x)+abs(behavior.events.map{tt}.y-y));

                    countMap{tt}(i,t,bb) = countMap{tt}(i,t,bb) + 1;
                    phaseMap{tt}{i} = [phaseMap{tt}{i}; bb t x y spikes.times{i}(sp(s))-start...
                        inst_rate(s) phases(round(spikes.times{i}(sp(s))*1250))];

                end
                end
            end
        end
        for m=1:length(trials)
            for mm = 1:length(behavior.events.trials{trials(m)}.mapping)
                [a b] = min(abs(behavior.events.trials{trials(m)}.x(mm)-behavior.events.map{tt}.x)+...
                    abs(behavior.events.trials{trials(m)}.y(mm)-behavior.events.map{tt}.y));
                occuMap{tt}(m,b) = occuMap{tt}(m,b) + 1;
            end
            occuMap{tt}(m,:) = medfilt1(occuMap{tt}(m,:),tau);
        end
    end
end

for tt =1:length(unique(behavior.events.trialConditions))
     trials = find(behavior.events.trialConditions==tt);
    for i = 1:length(spikes.times)
        numTrials = length(trials);
        for t = 1:length(trials)
%             rateMap{tt}(i,t,:) = smooth(squeeze(countMap{tt}(i,t,:))',tau)./ ...
%                 (smooth(occuMap{tt},tau)*(1/behavior.samplingRate)./numTrials);
            rateMap{tt}(i,t,:) = Smooth(squeeze(countMap{tt}(i,t,:))',tau)./ ...
                (Smooth(mean(occuMap{tt}),2*tau)*(1/behavior.samplingRate));
%               rateMap{tt}(i,t,:) = Smooth(squeeze(countMap{tt}(i,t,:))./...
%                   mean(occuMap{tt})'.*(1/behavior.samplingRate),tau);
        end
    end
end

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
firingMaps.sessionName = spikes.sessionName;
firingMaps.region = spikes.region; 

% backwards compatible but we'll eventually want to change this...
firingMaps.rateMaps = rateMap;
firingMaps.countMaps = countMap;
firingMaps.phaseMaps = phaseMap;
firingMaps.occupancy = occuMap;

if saveMat
   save([firingMaps.sessionName '.firingMaps.cellinfo.mat'],'firingMaps'); 
end



end
