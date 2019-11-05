function [firingMaps] = bz_firingMap1D(varargin)
% USAGE
% [firingMaps] = bz_firingMap1D(spikes,behavior,lfp,tau)
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%
%   behavior  - buzcode format behavior struct
%
%
%   tau       - desired smoothing window
%
%   saveMat   - logical (default: false) that saves firingMaps file
%
%
% OUTPUT
%
%   firingMaps - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%
% written by david tingley, 2017


% parse inputs
p=inputParser;
addRequired(p,'spikes',@isstruct);
addRequired(p,'behavior',@isstruct);
addRequired(p,'tau',@isnumeric);
addParameter(p,'saveMat',false,@islogical);
parse(p,varargin{:});

spikes = p.Results.spikes;
behavior = p.Results.behavior;
tau = p.Results.tau;
saveMat = p.Results.saveMat;

% start processing
for tt =1:length(unique(behavior.events.trialConditions))
    trials = find(behavior.events.trialConditions==tt);
    countMap{tt} = zeros(length(spikes.times),length(trials),length(behavior.events.map{tt}.x));
    occuMap{tt} = zeros(length(trials),length(behavior.events.map{tt}.x));
    if ~isempty(behavior.events.map{tt})

        for i =1:length(spikes.times)
            for t=1:length(trials)
                if ~isempty(spikes.times{i})
                start = behavior.events.trialIntervals(trials(t),1);
                stop =  behavior.events.trialIntervals(trials(t),2);
                f = find(spikes.times{i}>start);
                ff = find(spikes.times{i}<stop);
                sp = intersect(f,ff);
                for s =1:length(sp)
                    [a b] = min(abs(spikes.times{i}(sp(s))-behavior.events.trials{trials(t)}.timestamps(:,1)));

                    x = behavior.events.trials{trials(t)}.x(b);
                    y = behavior.events.trials{trials(t)}.y(b);
                    [aa bb] = min(abs(behavior.events.map{tt}.x-x)+abs(behavior.events.map{tt}.y-y));

                    countMap{tt}(i,t,bb) = countMap{tt}(i,t,bb) + 1;


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
            occuMap{tt}(m,:) = medfilt1(occuMap{tt}(m,:),5);
            occuMap{tt}(m,occuMap{tt}(m,:)<1) = 1;
        end
    end
end

for tt =1:length(unique(behavior.events.trialConditions))
    trials = find(behavior.events.trialConditions==tt);
    rateMap{tt} = zeros(length(spikes.times),length(trials),length(behavior.events.map{tt}.x));
    rateMap_box{tt} = zeros(length(spikes.times),length(trials),length(behavior.events.map{tt}.x));
    rateMap_unsmooth{tt} = zeros(length(spikes.times),length(trials),length(behavior.events.map{tt}.x));
    for i = 1:length(spikes.times)
        numTrials = length(trials);
        for t = 1:length(trials)
%             rateMap{tt}(i,t,:) = smooth(squeeze(countMap{tt}(i,t,:))',tau)./ ...
%                 (smooth(occuMap{tt},tau)*(1/behavior.samplingRate)./numTrials);
            rateMap{tt}(i,t,:) = Smooth(squeeze(countMap{tt}(i,t,:))',tau)./ ...
                (Smooth((occuMap{tt}(t,:)),tau)*(1/behavior.samplingRate));
            rateMap_box{tt}(i,t,:) = smooth(squeeze(countMap{tt}(i,t,:))',tau)./ ...
                (smooth((occuMap{tt}(t,:)),tau)*(1/behavior.samplingRate));
            rateMap_unsmooth{tt}(i,t,:) = (squeeze(countMap{tt}(i,t,:)))./ ...
                (Smooth((occuMap{tt}(t,:)),tau)*(1/behavior.samplingRate));
%               rateMap{tt}(i,t,:) = Smooth(squeeze(countMap{tt}(i,t,:))./...
%                   mean(occuMap{tt})'.*(1/behavior.samplingRate),tau);
        end
    end
    rateMap{tt}(isnan(rateMap{tt})) = 0;
    rateMap_box{tt}(isnan(rateMap_box{tt})) = 0;
    rateMap_unsmooth{tt}(isnan(rateMap_unsmooth{tt})) = 0;
end

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
firingMaps.sessionName = spikes.sessionName;
try
firingMaps.region = spikes.region; 
catch
   warning('spikes.region is missing') 
end
% backwards compatible but we'll eventually want to change this...
firingMaps.rateMaps = rateMap;
firingMaps.rateMaps_unsmooth = rateMap_unsmooth;
firingMaps.rateMaps_box = rateMap_box;
firingMaps.countMaps = countMap;
firingMaps.occupancy = occuMap;
firingMaps.dateRun = date;
if saveMat
   save([firingMaps.sessionName '.firingMaps.cellinfo.mat'],'firingMaps'); 
end

end
