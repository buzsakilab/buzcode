function [positionDecodingGLM] = bz_positionDecodingBayesian(varargin)
% USAGE
%   [] = bz_positionDecodingBayesian()
%
% INPUTS
%
%   spikes
%
%   behavior
%   
%   plotting   
%   
%   saveMat
%
% OUTPUTS
%
%   positionDecodingGLM     cellinfo struct with the following fields
%                      .UID
%                      .region
%                      .sessionName
%                      .results{1:nCells}
%                      .tau
%                      .condition
%
%
%
%
% SEE
%
% written by david tingley, 2017


%% parse inputs
p = inputParser();
addRequired(p,'spikes',@isstruct);
addRequired(p,'behavior',@isstruct);
addRequired(p,'lfp',@isstruct);
addRequired(p,'smoothingRange',@isvector);

addParameter(p,'plotting',@islogical);
addParameter(p,'saveMat',@islogical);
parse(p,varargin{:})

spikes = p.Results.spikes;
behavior = p.Results.behavior;
lfp = p.Results.lfp;
smoothingRange = p.Results.smoothingRange;
plotting = p.Results.plotting;
saveMat = p.Results.saveMat;


%% set up data format and cellinfo struct that will be returned

positionDecodingGLM.UID = spikes.UID;
positionDecodingGLM.region = spikes.region;
positionDecodingGLM.sessionName = spikes.sessionName; % session name
 


conditions = unique(behavior.events.trialConditions);
nCells = length(spikes.times);
positionSamplingRate = behavior.samplingRate;

% find a better way to get spike phase relationship...
[rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,5);

% iterate through conditions and compile spike trains and spike-phase
% trains
for cond = conditions
    trials = find(behavior.events.trialConditions==cond);
    intervals = behavior.events.trialIntervals(trials,:);
    
    for t = 1:length(trials)
        trial = trials(t);
        spk_trains{cond}{t} = zeros(nCells,ceil((intervals(t,2)-intervals(t,1))*1000)); % assumes intervals are in seconds, rounds to nearest millisecond
        phase_trains{cond}{t} = zeros(nCells,ceil((intervals(t,2)-intervals(t,1))*1000));
        for cell = 1:nCells
            if ~isempty(phaseMap{cond}{cell})
                f = find(phaseMap{cond}{cell}(:,2)==t);
                if ~isempty(f)
                for s=1:length(f)
                    phase_trains{cond}{t}(cell,ceil(phaseMap{cond}{cell}(f(s),5)*1000)) = ...
                        phaseMap{cond}{cell}(f(s),end);
                end
                end
            end
            sp = find(InIntervals(spikes.times{cell},intervals(t,:)));
            spk_trains{cond}{t}(cell,ceil((spikes.times{cell}(sp)-intervals(t,1))*1000))=1;
        end 
        
        position{cond}{t} = interp1(1:length(behavior.events.trials{trial}.x)...
            ,behavior.events.trials{trial}.mapping,1:positionSamplingRate/1000:length(...
            behavior.events.trials{trial}.x));
        position{cond}{t} = position{cond}{t}(1:length(spk_trains{cond}{t}));
    end
end

% collapse across trials..
for cond = conditions
   phase_trains{cond} = cell2mat(phase_trains{cond});
   spk_trains{cond} = cell2mat(spk_trains{cond});
   position{cond} = cell2mat(position{cond});
end


%% set up bayesian decoder



if saveMat 
   save([spikes.sessionName '.positionDecodingGLM.cellinfo.mat'],'positionDecodingGLM') 
end


