function [positionDecodingGLM] = bz_positionDecodingGLM(varargin)
% USAGE
%   [] = bz_positionDecodingGLM()
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

addParameter(p,'plotting',@islogical);
addParameter(p,'saveMat',@islogical);
parse(p,varargin{:})

spikes = p.Results.spikes;
behavior = p.Results.behavior;
lfp = p.Results.lfp;
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


%% set up GLM
for cell=1:nCells  % initialize table for data
    positionDecodingGLM.results{cell} = table;
end

for cond = conditions
    for window = smoothingRange
       for cell = 1:nCells 
            % phase coding 
            phase_trains_smooth=circ_smoothTS(phase_trains{cond}(cell,:),window,'method','mean','exclude',0); 
            [b dev stats] = glmfit(phase_trains_smooth',position{cond}','normal');
            yfit = glmval(b,phase_trains_smooth,'identity');
            struct.mse_phase = mean((yfit-position{cond}').^2); % mean squared error for rate code
            % phase coding 
            [b dev stats] = glmfit(cos(phase_trains_smooth)',position{cond}','normal');
            yfit = glmval(b,cos(phase_trains_smooth),'identity');
            struct.mse_phase_cos = mean((yfit-position{cond}').^2); % mean squared error for rate code
            % phase coding 
            [b dev stats] = glmfit(sin(phase_trains_smooth)',position{cond}','normal');
            yfit = glmval(b,sin(phase_trains_smooth),'identity');
            struct.mse_phase_sin = mean((yfit-position{cond}').^2); % mean squared error for rate code
            
            
            % rate coding
            rates_trains_smooth = smooth(spk_trains{cond}(cell,:),window);
            [b dev stats] = glmfit(rates_trains_smooth',position{cond}','normal');
            yfit = glmval(b,rates_trains_smooth,'identity');
            struct.mse_rate = mean((yfit-position{cond}').^2);  % mean squared error for rate code
            % extra variables to save
            struct.tau = window;
            struct.condition = cond;
            positionDecodingGLM.results{cell} = ...
                [positionDecodingGLM.results{cell};struct2table(struct)];
            if plotting & cell == 80
                subplot(2,2,1)
                plot(positionDecodingGLM.results{80}.tau,...
                    positionDecodingGLM.results{80}.mse_rate,'r')
                subplot(2,2,2)
                plot(positionDecodingGLM.results{80}.tau,...
                    positionDecodingGLM.results{80}.mse_phase_cos,'g')
                subplot(2,2,3)
                plot(positionDecodingGLM.results{80}.tau,...
                    positionDecodingGLM.results{80}.mse_phase_sin,'g')
                subplot(2,2,4)
                plot(positionDecodingGLM.results{80}.tau,...
                    positionDecodingGLM.results{80}.mse_phase,'g')
                pause(.0001)
            end
       end
    end
end


if saveMat 
   save([spikes.sessionName '.positionDecodingGLM.cellinfo.mat'],'positionDecodingGLM') 
end


