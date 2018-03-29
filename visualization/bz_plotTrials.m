function [] = bz_plotTrials(varargin)



p = inputParser;
addRequired(p,'behavior',@isstruct);
addParameter(p,'condition',[],@isnumeric);
addParameter(p,'color',[0 0 0],@isvector);
addParameter(p,'endPoints',true,@islogical);
addParameter(p,'spacing',1,@isnumeric);
parse(p,varargin{:})
behavior = p.Results.behavior;
condition = p.Results.condition;
color = p.Results.color;
endPoints = p.Results.endPoints;
spacing = p.Results.spacing;

if isempty(condition)
    f = factor(length(unique(behavior.events.trialConditions)));
    for tt = 1:length(unique(behavior.events.trialConditions))
        subplot(f(1),length(unique(behavior.events.trialConditions))./f(1),tt)
        ff = find(behavior.events.trialConditions==tt);

        for t = 1:length(ff)
        scatter(behavior.events.trials{ff(t)}.x(1:spacing:end),behavior.events.trials{ff(t)}.y(1:spacing:end),[],color,'.')
        hold on
        if endPoints
        scatter(behavior.events.trials{ff(t)}.x(1),behavior.events.trials{ff(t)}.y(1),'.g')
        scatter(behavior.events.trials{ff(t)}.x(end),behavior.events.trials{ff(t)}.y(end),'.r')
        end
        hold on
        end
    %     axis square
    title(['condition: ' num2str(tt) ', # of trials: ' num2str(length(ff))]);
    end
else
    for tt = 1:length(unique(behavior.events.trialConditions))
        ff = find(behavior.events.trialConditions==condition);

        for t = 1:length(ff)
        scatter(behavior.events.trials{ff(t)}.x(1:spacing:end),behavior.events.trials{ff(t)}.y(1:spacing:end),[],color,'.')
        hold on
        if endPoints
        scatter(behavior.events.trials{ff(t)}.x(1),behavior.events.trials{ff(t)}.y(1),'.g')
        scatter(behavior.events.trials{ff(t)}.x(end),behavior.events.trials{ff(t)}.y(end),'.r')
        hold on
        end
        end
    %     axis square
    title(['condition: ' num2str(condition) ', # of trials: ' num2str(length(ff))]);
    end
end

