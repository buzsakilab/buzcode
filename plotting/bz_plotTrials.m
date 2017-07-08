function [] = bz_plotTrials(behavior)



f = factor(length(unique(behavior.events.trialConditions)));

for tt = 1:length(unique(behavior.events.trialConditions))
    subplot(f(1),length(unique(behavior.events.trialConditions))./f(1),tt)
    ff = find(behavior.events.trialConditions==tt);
    
    for t = 1:length(ff)
        if t ~= 5 & t ~= 6
    scatter(behavior.events.trials{ff(t)}.x,behavior.events.trials{ff(t)}.y,'.k')
    hold on
    scatter(behavior.events.trials{ff(t)}.x(1),behavior.events.trials{ff(t)}.y(1),'.g')
    scatter(behavior.events.trials{ff(t)}.x(end),behavior.events.trials{ff(t)}.y(end),'.r')
    hold on
        end
    end
%     axis square
title(['condition: ' num2str(tt) ', # of trials: ' num2str(length(ff))]);
end

