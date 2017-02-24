function [sp_ale, sp] = SpindleSparsity(Sp, t_start, t_end )
% [sp_ale, sp] = SpindleSparsity(Sp, t_start, t_end )
% 
% computes the sparsity in spindle cycles
% INPUTS:
% Sp: a cell array of ts object (spike trains)
% t_start, t_end: the beginning and ends of the spindle cycles
%
% OUTPUTS:
% sp_ale: sparsity in the Treves definition \frac{< f  >^2}{<f^2>}
% sp: activity fraction (number of active cells), for each spindle cycle
  
%fpbatta 2003
  

  
  
  
  qs = zeros(length(Sp), length(t_start));
  
  
  for j = 1:length(Sp)
    h = Histogram_intervals(Data(Sp{j}), t_start, t_end);
    qs(j,:) = (h(:))';
  end
  
  
  f1 = sum(qs, 1) / length(Sp);
  f2 = sum(qs.*qs, 1) / length(Sp);
  
  sp_ale = (f1 .* f1) ./ f2;
  
  sp = sum((qs > 0), 1) / length(Sp);
  
  