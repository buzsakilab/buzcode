function [dm, ds] = DiagMean(M)
%dm = DiagMean(M) averages the diagonals of a matrix
%
% DiagMean computes the average of the diagonals of a matrix, usefule for
% example if you want to compute an average correlation as a function of
% time interval, from a time-bin correlation matrix
% INPUT:
% M: a square matrix
% OUTPUT:
% dm: the diagonal mean
% ds: standard deviation 
% batta 2001 status: under construction
  

  
  
if size(M, 1) ~= size(M, 2)
  error('M must be a square matrix');
end

l = size(M, 1);
dm = zeros(1, 2*l-1);

for i = -(l-1):(l-1)
  dm(l+i) = nanmean(diag(M, i));
  ds(l+i) = nanstd(diag(M,i));
end

  
  
  