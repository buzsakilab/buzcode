function [triplets_freq, triples_table] = spike_triplets(S, time_span)
%%  [triplets_freq, triples_table] = spike_triplets(S, time_span) 
%
% computes  a frequency table for all the spike triplets that occur
% within time_span
% INPUTS:
% S: a cell array of ts objects, containing spike trains
% time_span: the maximum duration spanned by a spike triplet
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
  
  
% mostly a wrapper for the count_triplets MEX function
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

[triplets_freq, triples_table] = count_triplets(t, idxs, n_cells, ...
						time_span);
