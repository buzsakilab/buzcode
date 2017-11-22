function [times groups spikesmat] = spikes2sorted(spikes)

times=[];
groups = [];
for i=1:length(spikes)
times = [times;spikes{i}];
groups = [groups;ones(length(spikes{i}),1)*i];
end
[a o] = sort(times);
times = times(o);
groups = groups(o);

spikesmat = [times,groups];

return