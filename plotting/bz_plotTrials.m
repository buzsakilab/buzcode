function [] = bz_plotTrials(trials)



f = factor(length(trials));

for tt = 1:length(trials)
subplot(length(trials)./f(1),f(1),tt)
for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
scatter(trials{tt}{t}(:,8),trials{tt}{t}(:,10),'.k')
hold on

%     axis([0 550 0 550])
end
for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
scatter(trials{tt}{t}(1,8),trials{tt}{t}(1,10),'.g')
scatter(trials{tt}{t}(end,8),trials{tt}{t}(end,10),'.r')
hold on

%     axis([0 550 0 550])
end
title(tt);
end

