function [phaseMaps] = bz_phaseMap1D(varargin)
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
                    
                    phaseMap{tt}{i} = [phaseMap{tt}{i}; bb t x y spikes.times{i}(sp(s))-start...
                        inst_rate(s) phases(round(spikes.times{i}(sp(s))*1250))];

                end
                end
            end
        end

    end
end

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
phaseMaps.UID = spikes.UID;
phaseMaps.sessionName = spikes.sessionName;
try
phaseMaps.region = spikes.region; 
catch
   warning('spikes.region is missing') 
end
% backwards compatible but we'll eventually want to change this...
phaseMaps.phaseMaps = phaseMap;
phaseMaps.dateRun = date;
if saveMat
   save([phaseMaps.sessionName '.phaseMaps.cellinfo.mat'],'phaseMaps'); 
end



end
