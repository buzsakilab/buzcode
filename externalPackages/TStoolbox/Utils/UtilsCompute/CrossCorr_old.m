function [C, B] = CrossCorr(t1, t2, binsize, nbins)
% [C, B] = CrossCorr(t1, t2, binsize, nbins)
%
% Cross Correlation of two time series
% INPUTS
% t1, t2: arrays containing sorted time series
% binsize: size of the bin for the cross correlation histogram
% nbins: number of bins in the histogram
% OUTPUTS
% C: the cross correlation histogram
% B: a vector with the time corresponding to the bins (optional)

% batta 1999
% MEX file
% status: beta