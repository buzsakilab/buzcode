function [idxs, range_min, range_max] = SplitQuantiles(D, n)
% [idxs, range_min, range_max] = SplitQuantiles(D, n) splits data in N
% quantiles 
%  
% INPUTS:
% D: a vector of data to be split 
% n: the desired number fo quantiles
% OUTPUTS:
% idxs: a cell array of index vectors, one per quantiles, from smaller to
% larger. Note that it's necessary that this be a cell array as it's not
% guaranteed that the quantiles will be of exactly the same size
% range_min, range_max: the minimum and maximum values in each quantile  
  
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  
  
  [D, ix] = sort(D);
  
  
  q_sz = floor(length(ix)/n);
  
  l = 0;
  
  for q = 1:n
    q_bottom = l+1;
    q_top = l+q_sz;
    l = q_top;
    if q == n
      q_top = length(ix);
    end
    
    
    range_min(q) = D(q_bottom);
    range_max(q) = D(q_top);
    idxs{q} = sort(ix(q_bottom:q_top));
  end
  
  
  
    