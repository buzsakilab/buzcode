function [positionDecodingBayesian] = bz_positionDecodingBayesian(varargin)
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
%   positionDecodingBayesian     cellinfo struct with the following fields
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

positionDecodingBayesian.UID = spikes.UID;
positionDecodingBayesian.region = spikes.region;
positionDecodingBayesian.sessionName = spikes.sessionName; % session name
 


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
            spk_trains{cond}{t}(cell,ceil((spikes.times{cell}(sp)-intervals(t,1))*1000+.00001))=1;
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

%% set up bayesian decoder
positionDecodingBayesian = table;

disp('running models...')
for cond = conditions
    for iter = 1:5
    r = randperm(length(phase_trains{cond}));
    for wind = smoothingRange
        for cell = 1:nCells
            phase_trains_smooth=[];
            rates_trains_smooth = [];
            position_train = [];
            for t = 1:length(phase_trains{cond})-1
                phase_trains_smooth=[phase_trains_smooth; circ_smoothTS(phase_trains{cond}{r(t)}(cell,:),wind,'method','mean','exclude',0)];
                rates_trains_smooth = [rates_trains_smooth; smooth(spk_trains{cond}{r(t)}(cell,:),wind)*wind];
                position_train = [position_train; position{cond}{r(t)}'];
            end
            phase_trains_smooth_train(cell,:) = phase_trains_smooth;
            rates_trains_smooth_train(cell,:) = rates_trains_smooth;
            phase_trains_smooth_test(cell,:)=[...
                circ_smoothTS(phase_trains{cond}{r(end)}(cell,:),wind,'method','mean','exclude',0)];
            rates_trains_smooth_test(cell,:) = [...
                smooth(spk_trains{cond}{r(end)}(cell,:),wind)*wind];
        end
        position_test = position{cond}{r(end)};
        
        %% rate coding model
        cl = poisson_naive_bayes_CL;
        cl = train(cl,round(rates_trains_smooth_train),position_train);
        
        for ts = 1:length(position_test)
           yfit_rate(ts) = test(cl,round(rates_trains_smooth_test(:,ts))); 
        end
        struct.mse_rate = mean((yfit_rate-position_test).^2);
        %% phase coding model
        
        % discretize phase_trains here...
        phase_trains_smooth_train_cos = cos(phase_trains_smooth_train);
        phase_trains_smooth_train_sin = sin(phase_trains_smooth_train);
        phase_trains_smooth_test_cos = cos(phase_trains_smooth_test);
        phase_trains_smooth_test_sin = sin(phase_trains_smooth_test);
        
        phase_trains_smooth_train(phase_trains_smooth_train==0)=nan;
        phase_trains_smooth_test(phase_trains_smooth_test==0)=nan;
        phase_trains_smooth_train = discretize(phase_trains_smooth_train,-pi:.1:pi);
        phase_trains_smooth_test = discretize(phase_trains_smooth_test,-pi:.1:pi);
        phase_trains_smooth_train(isnan(phase_trains_smooth_train))=0;
        phase_trains_smooth_test(isnan(phase_trains_smooth_test))=0;
        
        phase_trains_smooth_train_cos(phase_trains_smooth_train_cos==0)=nan;
        phase_trains_smooth_test_cos(phase_trains_smooth_test_cos==0)=nan;
        phase_trains_smooth_train_cos = discretize(phase_trains_smooth_train_cos,-pi:.1:pi);
        phase_trains_smooth_test_cos = discretize(phase_trains_smooth_test_cos,-pi:.1:pi);
        phase_trains_smooth_train_cos(isnan(phase_trains_smooth_train_cos))=0;
        phase_trains_smooth_test_cos(isnan(phase_trains_smooth_test_cos))=0;
        
        phase_trains_smooth_train_sin(phase_trains_smooth_train_sin==0)=nan;
        phase_trains_smooth_test_sin(phase_trains_smooth_test_sin==0)=nan;
        phase_trains_smooth_train_sin = discretize(phase_trains_smooth_train_sin,-pi:.1:pi);
        phase_trains_smooth_test_sin = discretize(phase_trains_smooth_test_sin,-pi:.1:pi);
        phase_trains_smooth_train_sin(isnan(phase_trains_smooth_train_sin))=0;
        phase_trains_smooth_test_sin(isnan(phase_trains_smooth_test_sin))=0;
        
        % non-transformed circular decoding
        cl = poisson_naive_bayes_CL;
        cl = train(cl,phase_trains_smooth_train,position_train);
        
        for ts = 1:length(position_test)
           yfit_circ(ts) = test(cl,phase_trains_smooth_test(:,ts)); 
        end
        struct.mse_phase = mean((yfit_circ-position_test).^2);
        
        % cos and sin transformed circular decoding
        cl = poisson_naive_bayes_CL;
        cl = train(cl,phase_trains_smooth_train_cos,position_train);
        
        for ts = 1:length(position_test)
           yfit_cos(ts) = test(cl,phase_trains_smooth_test_cos(:,ts)); 
        end
        struct.mse_phase_cos = mean((yfit_cos-position_test).^2);
        
        cl = poisson_naive_bayes_CL;
        cl = train(cl,phase_trains_smooth_train_sin,position_train);
        
        for ts = 1:length(position_test)
           yfit_sin(ts) = test(cl,phase_trains_smooth_test_sin(:,ts)); 
        end
        struct.mse_phase_sin = mean((yfit_sin-position_test).^2);
        
        
        %% put data into struct/table
        struct.tau = wind;
        struct.condition = cond;
        struct.iter = iter;
        struct.trialOrder = r;
        positionDecodingBayesian = [positionDecodingBayesian;struct2table(struct)];
        if plotting
            subplot(2,2,1)
            plot(positionDecodingBayesian.tau,...
                positionDecodingBayesian.mse_rate,'r')
            subplot(2,2,2)
            plot(positionDecodingBayesian.tau,...
                positionDecodingBayesian.mse_phase_cos,'g')
            subplot(2,2,3)
            plot(positionDecodingBayesian.tau,...
                positionDecodingBayesian.mse_phase_sin,'g')
            subplot(2,2,4)
            plot(positionDecodingBayesian.tau,...
                positionDecodingBayesian.mse_phase,'g')
            pause(.001)
        end
        clear *train *test
        disp(['finished with window: ' num2str(wind) ' out of ' num2str(smoothingRange(end)) ' total'])
    end
    end
    disp(['finished with condition: ' num2str(cond) ' out of ' num2str(length(conditions)) ' total'])
end


positionDecodingBayesion.dataRun = date;

if saveMat 
   save([spikes.sessionName '.positionDecodingBayesian.cellinfo.mat'],'positionDecodingBayesian') 
end


