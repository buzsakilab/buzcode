function [positionDecodingGLM_gaussian] = bz_positionDecodingGLM_gauss(varargin)
% USAGE
%   [] = bz_positionDecodingGLM_gauss()
%
% INPUTS
%
%   spikes
%
%   behavior
%
%   lfp 
%
%   smoothingRange
%   
%   plotting   - cell UID to plot
%   
%   saveMat
%
% OUTPUTS
%
%   positionDecodingGLM_gaussian     cellinfo struct with the following fields
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

addParameter(p,'plotting',true,@isnumeric);
addParameter(p,'saveMat',[],@islogical);
parse(p,varargin{:})

spikes = p.Results.spikes;
behavior = p.Results.behavior;
lfp = p.Results.lfp;
smoothingRange = p.Results.smoothingRange;
plotting = p.Results.plotting;
saveMat = p.Results.saveMat;


%% set up data format and cellinfo struct that will be returned
disp('initializing data structs and formatting input data...')

positionDecodingGLM_gaussian.UID = spikes.UID;
positionDecodingGLM_gaussian.region = spikes.region;
positionDecodingGLM_gaussian.sessionName = spikes.sessionName; % session name
 
conditions = unique(behavior.events.trialConditions);
nCells = length(spikes.times);
positionSamplingRate = behavior.samplingRate;

% find a better way to get spike phase relationship...
[firingMaps] = bz_firingMap1D(spikes,behavior,lfp,5);

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
            if ~isempty(firingMaps.phaseMaps{cond}{cell})
                f = find(firingMaps.phaseMaps{cond}{cell}(:,2)==t);
                if ~isempty(f)
                for s=1:length(f)
                    phase_trains{cond}{t}(cell,ceil(firingMaps.phaseMaps{cond}{cell}(f(s),5)*1000)) = ...
                        firingMaps.phaseMaps{cond}{cell}(f(s),end);
                end
                end
            end
            sp = find(InIntervals(spikes.times{cell},intervals(t,:)));
            if ~isempty(sp)
                spks = Restrict(ceil((spikes.times{cell}(sp)-intervals(t,1))*1000+.0000001),[1 size(spk_trains{cond}{t},2)]);
                spk_trains{cond}{t}(cell,spks)=1;
            end
        end 
        
%         position{cond}{t} = interp1(1:length(behavior.events.trials{trial}.x)...
%             ,behavior.events.trials{trial}.mapping,1:positionSamplingRate/1000:length(...
%             behavior.events.trials{trial}.x));
        nBins = size(spk_trains{cond}{t},2);
        nPos = length(behavior.events.trials{trial}.x);
        position{cond}{t} = interp1(1:length(behavior.events.trials{trial}.x)...
            ,behavior.events.trials{trial}.mapping,1:(nPos-1)/nBins:length(...
            behavior.events.trials{trial}.x)); % the -1 gaurantees the length to be longer than the above spk/phase trains
        position{cond}{t} = position{cond}{t}(1:length(spk_trains{cond}{t}));
        
    end
end

% collapse across trials..
for cond = conditions
%    phase_trains{cond} = cell2mat(phase_trains{cond});
%    spk_trains{cond} = cell2mat(spk_trains{cond});
% this was a bad idea because it lead to smoothing across trial boundaries,
% fixed below by smoothing first, then concatenating trials...
   position{cond} = cell2mat(position{cond});
end


%% set up GLM
for cell=1:nCells  % initialize table for data
    positionDecodingGLM_gaussian.results{cell} = table;
end
disp('running models...')
for cond = conditions
    warning off
    for wind = smoothingRange
       for cell = 80%1:nCells 
            %smoothing
            phase_trains_smooth=[];
            cos_phase_trains_smooth=[];
            sin_phase_trains_smooth=[];
            rates_trains_smooth = [];
            for t = 1:length(phase_trains{cond})
%                 phase_trains_smooth=[phase_trains_smooth;...
%                     circ_smoothTS(phase_trains{cond}{t}(cell,:),wind,'method','mean','exclude',0)];
                phase_trains_smooth=[phase_trains_smooth;...
                     gaussfilt(1:length(phase_trains{cond}{t}(cell,:)),phase_trains{cond}{t}(cell,:),wind)'];
                cos_phase_trains_smooth=[cos_phase_trains_smooth;...
                     gaussfilt(1:length(phase_trains{cond}{t}(cell,:)),cos(phase_trains{cond}{t}(cell,:)),wind)'];
                sin_phase_trains_smooth=[sin_phase_trains_smooth;...
                     gaussfilt(1:length(phase_trains{cond}{t}(cell,:)),sin(phase_trains{cond}{t}(cell,:)),wind)'];
                rates_trains_smooth = [rates_trains_smooth; ...
                                       gaussfilt(1:length(phase_trains{cond}{t}(cell,:)),spk_trains{cond}{t}(cell,:),wind)'];
            end
            % phase coding
            [b dev stats] = glmfit(phase_trains_smooth',position{cond}','normal');
            yfit = glmval(b,phase_trains_smooth,'identity');
            struct.mse_phase = mean((yfit-position{cond}').^2); % mean squared error for rate code
            struct.mse_phase_pval = stats.p';

            % phase coding 
            [b dev stats] = glmfit(cos_phase_trains_smooth',position{cond}','normal');
            yfit = glmval(b,cos_phase_trains_smooth,'identity');
            struct.mse_phase_cos = mean((yfit-position{cond}').^2); % mean squared error for rate code
            struct.mse_phase_cos_pval = stats.p';

            % phase coding 
            [b dev stats] = glmfit(sin_phase_trains_smooth',position{cond}','normal');
            yfit = glmval(b,sin_phase_trains_smooth,'identity');
            struct.mse_phase_sin = mean((yfit-position{cond}').^2); % mean squared error for rate code
            struct.mse_phase_sin_pval = stats.p';
            
            % phase coding all
            [b dev stats] = glmfit([sin_phase_trains_smooth,cos_phase_trains_smooth],position{cond}','normal');
            yfit = glmval(b,[sin_phase_trains_smooth,cos_phase_trains_smooth],'identity');
            struct.mse_phase_all = mean((yfit-position{cond}').^2); % mean squared error for rate code
            struct.mse_phase_all_pval = stats.p';
            
            % rate coding
            [b dev stats] = glmfit(rates_trains_smooth',position{cond}','normal');
            yfit = glmval(b,rates_trains_smooth,'identity');
            struct.mse_rate = mean((yfit-position{cond}').^2);  % mean squared error for rate code
            struct.mse_rate_pval = stats.p';
            
            % extra variables to save
            struct.dfe = stats.dfe;
            struct.tau = wind;
            struct.condition = cond;
            positionDecodingGLM_gaussian.results{cell} = ...
                [positionDecodingGLM_gaussian.results{cell};struct2table(struct)];
            
            if ~isempty(plotting) & cell == plotting  % can only display one neuron for now..
%                 figure(cond);
                subplot(2,2,1)
                title('GLM decoding of pos, r-rate, g-phase')
                rows = positionDecodingGLM_gaussian.results{plotting}.condition == cond;
                plot(positionDecodingGLM_gaussian.results{plotting}.tau(rows),...
                    positionDecodingGLM_gaussian.results{plotting}.mse_rate(rows),'r')
                hold on
                plot(positionDecodingGLM_gaussian.results{plotting}.tau(rows),...
                    positionDecodingGLM_gaussian.results{plotting}.mse_phase_all(rows),'g')
                hold off
                subplot(2,2,2)
                plot(positionDecodingGLM_gaussian.results{plotting}.tau(rows),...
                    positionDecodingGLM_gaussian.results{plotting}.mse_phase_cos(rows),'g')
                subplot(2,2,3)
                plot(positionDecodingGLM_gaussian.results{plotting}.tau(rows),...
                    positionDecodingGLM_gaussian.results{plotting}.mse_phase_sin(rows),'g')
                subplot(2,2,4)
                plot(positionDecodingGLM_gaussian.results{plotting}.tau(rows),...
                    positionDecodingGLM_gaussian.results{plotting}.mse_phase(rows),'g')
                pause(.0001)
            end
       end
       disp(['finished with wind: ' num2str(wind) ' out of ' num2str(smoothingRange(end)) ' total'])
    end
    positionDecodingGLM_gaussian.dateRun = date;  % this can take a very long time so lets save each loop...
    if saveMat 
        save([spikes.sessionName '.positionDecodingGLM_gaussian.cellinfo.mat'],'positionDecodingGLM_gaussian') 
    end
    disp(['finished with condition: ' num2str(cond) ' out of ' num2str(length(conditions)) ' total'])
end






