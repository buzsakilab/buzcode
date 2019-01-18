function [times groups spikesmat] = spikes2sorted(spikes)

times=[];
groups = [];
for i=1:length(spikes)
%     if ~isempty(spikes{i})
        times = [times;spikes{i}];
        groups = [groups;ones(length(spikes{i}),1)*i];
%     else
%         times = [times; nan];
%         groups = [groups; i];
%     end
end

[a o] = sort(times);
times = times(o);
groups = groups(o);

spikesmat = [times,groups];

return