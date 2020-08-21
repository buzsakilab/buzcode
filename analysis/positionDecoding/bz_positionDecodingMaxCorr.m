function [positionDecodingMaxCorr] = bz_positionDecodingMaxCorr(varargin)
% USAGE
%   [] = bz_positionDecodingMaxCorr()
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
%   plotting   
%   
%   saveMat
%
% OUTPUTS
%
%   positionDecodingMaxCorr     cellinfo struct with the following fields
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


%% set up data format and popinfo struct that will be returned


positionDecodingMaxCorr.UID = spikes.UID;
positionDecodingMaxCorr.region = spikes.region;
positionDecodingMaxCorr.sessionName = spikes.sessionName; % session name

conditions = unique(behavior.events.trialConditions);
nCells = length(spikes.times);
positionSamplingRate = behavior.samplingRate;

% find a better way to get spike phase relationship...
[firingMaps] = bz_firingMap1D(spikes,behavior,5);
[phaseMaps] = bz_phaseMap1D(spikes,behavior,lfp,5);
[b a] = butter(3,[4/625 12/625],'bandpass');
theta_phases = angle(hilbert(FiltFiltM(b,a,double(lfp.data))));
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
            if ~isempty(phaseMaps.phaseMaps{cond}{cell})
                f = find(phaseMaps.phaseMaps{cond}{cell}(:,2)==t);
                if ~isempty(f)
                for s=1:length(f)
                    phase_trains{cond}{t}(cell,ceil(phaseMaps.phaseMaps{cond}{cell}(f(s),5)*1000)) = ...
                        phaseMaps.phaseMaps{cond}{cell}(f(s),end);
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
        take = FindInInterval(lfp.timestamps,intervals(t,:));
        theta{cond}{t} = makeLength(theta_phases(take(1):take(2)),length(position{cond}{t}));
    end
end

%% set up MaxCorr decoder
positionDecodingMaxCorr.results = table;

disp('running models...')
warning off % max cl model returns warning with NaNs and Inf values..
for cond = conditions
    r = randperm(length(phase_trains{cond}));
    for iter = 1:5%length(r)
        r = circshift(r,1);
    for wind = smoothingRange
        for cell = 1:nCells
            phase_trains_smooth=[];
            rates_trains_smooth = [];
            position_train = [];
            theta_train = [];
            for t = 1:length(phase_trains{cond})-1
                phase_trains_smooth=[phase_trains_smooth; circ_smoothTS(phase_trains{cond}{r(t)}(cell,:),wind,'method','mean','exclude',0)];
                rates_trains_smooth = [rates_trains_smooth; smooth(spk_trains{cond}{r(t)}(cell,:),wind)*wind];
                position_train = [position_train; position{cond}{r(t)}'];
                theta_train = [theta_train, theta{cond}{r(t)}];
            end
            phase_trains_smooth_train(cell,:) = phase_trains_smooth;
            rates_trains_smooth_train(cell,:) = rates_trains_smooth;
            phase_trains_smooth_test(cell,:)=[...
                circ_smoothTS(phase_trains{cond}{r(end)}(cell,:),wind,'method','mean','exclude',0)];
            rates_trains_smooth_test(cell,:) = [...
                smooth(spk_trains{cond}{r(end)}(cell,:),wind)*wind];
            
        end
        position_test = position{cond}{r(end)};
        theta_test = theta{cond}{r(end)};
        
        %% rate coding model
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[rates_trains_smooth_train; theta_train],round(position_train));
        
        for ts = 1:length(position_test)
           yfit_rate(ts) = test(cl,[round(rates_trains_smooth_test(:,ts));theta_test(ts)]); 
        end
        struct.mse_rate = mean((yfit_rate-position_test).^2);
        % chance rate
        rr = randperm(length(theta_train));
        rrr = randperm(length(theta_test));
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[rates_trains_smooth_train(:,rr); theta_train],round(position_train));
        
        for ts = 1:length(position_test)
           yfit_chance_rate(ts) = test(cl,[round(rates_trains_smooth_test(:,rrr(ts)));theta_test(ts)]); 
        end
        struct.mse_chance_rate = mean((yfit_chance_rate-position_test).^2);
        
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
        phase_trains_smooth_train_cos = discretize(phase_trains_smooth_train_cos,-1:.1:1);
        phase_trains_smooth_test_cos = discretize(phase_trains_smooth_test_cos,-1:.1:1);
        phase_trains_smooth_train_cos(isnan(phase_trains_smooth_train_cos))=0;
        phase_trains_smooth_test_cos(isnan(phase_trains_smooth_test_cos))=0;
        
        phase_trains_smooth_train_sin(phase_trains_smooth_train_sin==0)=nan;
        phase_trains_smooth_test_sin(phase_trains_smooth_test_sin==0)=nan;
        phase_trains_smooth_train_sin = discretize(phase_trains_smooth_train_sin,-1:.1:1);
        phase_trains_smooth_test_sin = discretize(phase_trains_smooth_test_sin,-1:.1:1);
        phase_trains_smooth_train_sin(isnan(phase_trains_smooth_train_sin))=0;
        phase_trains_smooth_test_sin(isnan(phase_trains_smooth_test_sin))=0;
        
        % non-transformed circular decoding
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[phase_trains_smooth_train; theta_train],round(position_train));
        
        for ts = 1:length(position_test)
           yfit_circ(ts) = test(cl,[phase_trains_smooth_test(:,ts);theta_test(ts)]); 
        end
        struct.mse_phase = mean((yfit_circ-position_test).^2);
        
        % cos and sin transformed circular decoding
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[phase_trains_smooth_train_cos; theta_train],round(position_train));
        
        for ts = 1:length(position_test)
           yfit_cos(ts) = test(cl,[phase_trains_smooth_test_cos(:,ts);theta_test(ts)]); 
        end
        struct.mse_phase_cos = mean((yfit_cos-position_test).^2);
        
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[phase_trains_smooth_train_sin; theta_train],round(position_train));
        
        for ts = 1:length(position_test)
           yfit_sin(ts) = test(cl,[phase_trains_smooth_test_sin(:,ts);theta_test(ts)]); 
        end
        struct.mse_phase_sin = mean((yfit_sin-position_test).^2);
        
        % all phase models in one...
        cl = max_correlation_coefficient_CL;
        all_phase_train = [phase_trains_smooth_train_cos;phase_trains_smooth_train_sin];
        all_phase_test =  [phase_trains_smooth_test_cos;phase_trains_smooth_test_sin];
        cl = train(cl,[all_phase_train; theta_train],round(position_train));
        
        for ts = 1:length(position_test)
           yfit_circ_all(ts) = test(cl,[all_phase_test(:,ts);theta_test(ts)]); 
        end
        struct.mse_phase_all = mean((yfit_circ_all-position_test).^2);
        
        % chance phase
        cl = max_correlation_coefficient_CL;
        all_phase_train = [phase_trains_smooth_train_cos;phase_trains_smooth_train_sin];
        all_phase_test =  [phase_trains_smooth_test_cos;phase_trains_smooth_test_sin];
        cl = train(cl,[all_phase_train(:,rr); theta_train],round(position_train));
        
        for ts = 1:length(position_test)
           yfit_chance(ts) = test(cl,[all_phase_test(:,rrr(ts));theta_test(ts)]); 
        end
        struct.mse_chance = mean((yfit_chance-position_test).^2);
        
        %% put data into struct/table
        struct.tau = wind;
        struct.condition = cond;
        struct.iter = iter;
        struct.trialOrder = r;
        positionDecodingMaxCorr.results = [positionDecodingMaxCorr.results;struct2table(struct)];
        if plotting
            clf
            subplot(2,2,1)
            t_rate = varfun(@mean,positionDecodingMaxCorr.results,'InputVariables','mse_rate',...
            'GroupingVariables',{'tau','condition'});
            t_phase = varfun(@mean,positionDecodingMaxCorr.results,'InputVariables','mse_phase',...
            'GroupingVariables',{'tau','condition'});
            t_phase_cos = varfun(@mean,positionDecodingMaxCorr.results,'InputVariables','mse_phase_cos',...
            'GroupingVariables',{'tau','condition'});
            t_phase_sin = varfun(@mean,positionDecodingMaxCorr.results,'InputVariables','mse_phase_sin',...
            'GroupingVariables',{'tau','condition'});
            t_chance = varfun(@mean,positionDecodingMaxCorr.results,'InputVariables','mse_chance',...
            'GroupingVariables',{'tau','condition'});
            t_chance_s = varfun(@std,positionDecodingMaxCorr.results,'InputVariables','mse_chance',...
            'GroupingVariables',{'tau','condition'});
            t_chance_rate = varfun(@mean,positionDecodingMaxCorr.results,'InputVariables','mse_chance_rate',...
            'GroupingVariables',{'tau','condition'});
            t_chance_s_rate = varfun(@std,positionDecodingMaxCorr.results,'InputVariables','mse_chance_rate',...
            'GroupingVariables',{'tau','condition'});
            t_rate_s = varfun(@std,positionDecodingMaxCorr.results,'InputVariables','mse_rate',...
            'GroupingVariables',{'tau','condition'});
            t_phase_s = varfun(@std,positionDecodingMaxCorr.results,'InputVariables','mse_phase',...
            'GroupingVariables',{'tau','condition'});
            t_phase_cos_s = varfun(@std,positionDecodingMaxCorr.results,'InputVariables','mse_phase_cos',...
            'GroupingVariables',{'tau','condition'});
            t_phase_sin_s = varfun(@std,positionDecodingMaxCorr.results,'InputVariables','mse_phase_sin',...
            'GroupingVariables',{'tau','condition'});
            tab = join(join(join(join(t_rate,t_phase),t_chance),t_phase_cos),t_phase_sin);
            tab_s = join(join(join(join(t_rate_s,t_phase_s),t_chance_s),t_phase_sin_s),t_phase_cos_s);
            rows = find(tab.condition==cond);

            title('MaxCorr decoding of pos, r-rate, g-phase')
%             iterRows = positionDecodingMaxCorr.iter == iter;
%             plot(positionDecodingMaxCorr.results.tau,...
%                 positionDecodingMaxCorr.results.mse_rate,'r')
%             hold on
%             plot(positionDecodingMaxCorr.results.tau,...
%                 positionDecodingMaxCorr.results.mse_phase_all,'g')
%             hold off
%             boundedline(tab.tau,t_chance_rate.mean_mse_chance_rate(rows),t_chance_s_rate.std_mse_chance_rate(rows),'k')
            boundedline(tab.tau,tab.mean_mse_chance(rows),tab_s.std_mse_chance(rows),'k')
            boundedline(tab.tau,tab.mean_mse_phase_cos(rows),tab_s.std_mse_phase_cos(rows),'g')
            boundedline(tab.tau,tab.mean_mse_rate(rows),tab_s.std_mse_rate(rows),'r')
%             set(gca,'xscale','log')
            subplot(2,2,2)
            boundedline(tab.tau,tab.mean_mse_phase(rows),tab_s.std_mse_phase(rows),'.g')
            boundedline(tab.tau,tab.mean_mse_phase_sin(rows),tab_s.std_mse_phase_sin(rows),'--g')
            boundedline(tab.tau,tab.mean_mse_phase_cos(rows),tab_s.std_mse_phase_cos(rows),'g')
%             subplot(2,2,3)
%             plot(positionDecodingMaxCorr.results.tau,...
%                 positionDecodingMaxCorr.results.mse_phase_sin,'g')
%             subplot(2,2,4)
%             plot(positionDecodingMaxCorr.results.tau,...
%                 positionDecodingMaxCorr.results.mse_phase,'g')
            pause(.001)
        end
        clear *train *test yfit*
        disp(['finished with window: ' num2str(wind) ' out of ' num2str(smoothingRange(end)) ' total'])
    end
    end
    disp(['finished with condition: ' num2str(cond) ' out of ' num2str(length(conditions)) ' total'])
    positionDecodingMaxCorr.dateRun = date;
    if saveMat 
       save([spikes.sessionName '.positionDecodingMaxCorr.popinfo.mat'],'positionDecodingMaxCorr') 
    end
end





