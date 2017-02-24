function rank = RankSpikes(S, threshold)
% rank = RankSpikes(S, threshold)
%
% gives each spike a "number" representing the number from the entrance in
% the place field
%
% INPUTS:
% S: a cell array of ts spike trains.
% the time threshold to decide where to set the boundaries between spikes
% relative to different place field visits.
% OUTPUT:
% rank: cell array of tsd objects

% batta 2000
% status: alpha


for iC = 1:length(S)
  s = Data(S{iC});
  sd = diff(s);
  rk = zeros(size(S{iC}));
  rk(1) = 1;
  for i = 2:length(s)
    if sd(i-1) > threshold
      rk(i) = 1;
    else
      rk(i) = rk(i-1)+1;
    end
  end
  rank{iC} = tsd(s,rk);
end
