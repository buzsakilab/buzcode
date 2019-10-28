function [order, excluded] = bz_plotRasterTrial(spikes,interval,countLimit,order,ex)


for c= 1:length(spikes.times)
    [idx] =find(InIntervals(spikes.times{c},interval));
    ts{c} = spikes.times{c}(idx);
    mu(c) = median(ts{c});
end

% mu(isnan(mu)) = 0;
if isempty(order)
    [a b] = sort(mu');
    order = b;
else
    b = order;
end

excluded = zeros(length(spikes.times),1);

for c= 1:length(spikes.times)
    if length(ts{b(c)}) < countLimit & ex(c) ~= 1
%         scatter(ts{b(c)}-interval(1),ones(length(ts{b(c)}),1)*c,'.')
        scatter(ts{b(c)},ones(length(ts{b(c)}),1)*c,'.')
        hold on
    else
        scatter(interval(1)-.001,c,'.')
        hold on
        excluded(c) = 1;
    end
end

xlim([interval])
ylabel('cell #')
xlabel('time (s)')