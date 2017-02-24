function S = SpkTrainGenerator(r,dt)

% S = SpkTrainGenerator(r,time)

r = r(:);

rng('shuffle') % Seed random number generator

% We want a mean firing rate of r Hz. So we want r spikes in
% 1 second. There are 1/dt points in 1 second so we want r
% spikes every 1/dt points.
pspike = r*dt; % Probability of Spiking at each point in time
S = pspike > rand(size(r));
