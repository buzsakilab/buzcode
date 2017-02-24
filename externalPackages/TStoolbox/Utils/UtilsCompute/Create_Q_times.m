function Q = Create_Q_times(data_set, tsmatrix, dt)
% Q = Create_Q_times(data_set, tsmatrix, dt)
% Purpose: Create a sparse Q matrix
%
% Input: a cell array of ts objects containing the spike trains
% start and end timestamps, and a dt interval.
%
% Output: a tsd object containing a sparse Q matrix  as data
%


if nargin == 2, dt = 100, end; % bin size in milliseconds

% Load in the data sets.


S = data_set;
Q_tsarray = MakeQfromS(S, dt*10);
Q_tsarray = tsd(Range(Q_tsarray, 'ts'), Data(Q_tsarray));
% Find Alignment finds the indices in the Q matrix of the timestamps.
idx_i = findAlignment(Q_tsarray, tsmatrix(:,1));
idx_j = findAlignment(Q_tsarray, tsmatrix(:,2));
% Get rid of the negative time stamps. Why are they negative? I think
% because they are not aligned.
if length(find(idx_i<0)) > 0 | length(find(idx_j<0)) > 0
  disp('ERROR, dumping indices..');
  [idx_i idx_j]
  error('There were negative index values in your timstamp matrix.');
  error('Your interval periods are wrong or your tstamps are bad.');
end
disp('Extracting intervals');
Q = Extract_intervals_times(Q_tsarray, idx_i, idx_j);
