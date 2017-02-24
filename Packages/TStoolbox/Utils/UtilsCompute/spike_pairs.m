function [pairs_freq, pairs_table, pairs_ts] = ...
    spike_pairs(S,  time_span, on_same_tet)
%%  [pairs_freq, pairs_table] = spike_pairs(S, time_span) 
%
% computes  a frequency table for all the spike pairs that occur
% within time_span
% INPUTS:
% S: a cell aray of ts objects, containing spike trains
% time_span: the maximum duration spanned by a spike triplet
% on_same_tet: matrix indicating cells from same trode
% OUTPUTS:
% pairs_freq: a n_comb(N_cells, 2) X 2 matrix of frequencies, rows a
% are cell pairs columns possible orderings, according to 
% 1: 1 2
% 2: 2 1
% pairs_table: a n_comb(N_cells, 2) X 2 matrix, a lookup table for the
% row dimension in pairs_freq, i-th row is the set of cell indices for
% the i-th triple, in increasing order
% pairs_ts
%
% mostly a wrapper for the count_pairs_ts MEX function
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




[pairs_freq, pairs_table, pairs_ts] = ...
    count_pairs(t, idxs, n_cells, ...
		      time_span, on_same_tet);
