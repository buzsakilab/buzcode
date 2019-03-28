function bz_plotRasterTrial(spikes,interval,countLimit)


for c= 1:length(spikes.times)
[idx] =find(InIntervals(spikes.times{c},interval));
ts{c} = spikes.times{c}(idx);
mu(c) = median(ts{c});
end

% mu(isnan(mu)) = 0;
[a b] = sort(mu');

for c= 1:length(spikes.times)
    if length(ts{b(c)}) < countLimit
    scatter(ts{b(c)}-interval(1),ones(length(ts{b(c)}),1)*c,'.')
    hold on
    end
end

ylabel('cell #')
xlabel('time (s)')