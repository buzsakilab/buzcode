function [triplets_freq, triples_table, triplets_ts, triplets_x, triplets_y, triplets_phi] = ...
    spike_triplets_ts(S, XS, YS, bins_maze, time_span, on_same_tet)
%%  [triplets_freq, triples_table] = spike_triplets(S, time_span) 
%
% computes  a frequency table for all the spike triplets that occur
% within time_span
% INPUTS:
% S: a cell aray of ts objects, containing spike trains
% XS
% YS
% bins_maze: 3 tsd objs containing x, y, and angular coords
% time_span: the maximum duration spanned by a spike triplet
% on_same_tet: matrix indicating cells from same trode
% OUTPUTS:
% triplets_freq: a n_comb(N_cells, 3) X 6 matrix of frequencies, rows a
% are cell triplets columns possible orderings, according to 
% 1: 1 2 3
% 2: 1 3 2
% 3: 2 1 3 
% 4: 2 3 1
% 5: 3 1 2 
% 6: 3 2 1
% triples_table: a n_comb(N_cells, 3) X 3 matrix, a lookup table for the
% row dimension in triplets_freq, i-th row is the set of cell indices for
% the i-th triple, in increasing order
% triplets_ts
% triplets_x
% triplets_y
% triplets_phi: 4 n_comb(N_cells, 3) X 6 cell array of vectors containing
% the timestamps, x, y, and angular coords at which each triplet occurred 

% mostly a wrapper for the count_triplets_ts MEX function
% batta 2002, in collaboration with A. Treves
  
  
  

n_cells = length(S);


t = zeros(0,1);
idxs = zeros(0,1);
for i = 1:length(S)
  t = [t ; Data(S{i})];
  idxs = [idxs; i * ones(size(Data(S{i})))];
end


[t, ix] = sort(t);
idxs = idxs(ix);
x = Restrict(XS, t);
y = Restrict(YS, t);
phi = Restrict(bins_maze, t);

p = Data(phi);
p(find(t < StartTime(bins_maze) | t > EndTime(bins_maze))) = -1;
phi = tsd(t, p);



[triplets_freq, triples_table, triplets_ts, triplets_x, triplets_y, triplets_phi] = ...
    count_triplets_ts(t, idxs, Data(x), Data(y), Data(phi), n_cells, ...
		      time_span, on_same_tet);
