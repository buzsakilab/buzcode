function SE = CR2SE(cr, threshold)

sz = size(cr);

if(sz(1) > 1)
  cr = cr';
end

l = length(cr);

d1 = [0 diff(cr)];
d = [diff(cr) 0];
cr1 = [cr 0];
cr2 = [0 cr];

cross_point = find(cr1 > threshold & cr2 < threshold);

maxes = find(d1 > 0 & d < 0);


spike_peak = maxes(FindAlignNext(maxes, cross_point));

spike_peak = spike_peak';
nspikes = length(cross_point);

sp_point = repmat(spike_peak, 1, 32) + repmat(-7:24, nspikes, 1);
%sp_point = repmat(cross_point', 1, 32) + repmat(-7:24, nspikes, 1);

spikes = cr(sp_point);
SE = tsd(spike_peak, spikes);

