%  [C, B] = CrossCorr(t1,t2,binsize,nbins);
%
%  cross correlations
%  MEX file
% 
% batta 1999
% 
% input: t1, t2: two time series to cross correlate in timestamp units (1/10000 sec)
% (assumed to be sorted) 
%        binsize: the binsize for the cross corr histogram in msec
%        nbins: the number of bins
% output: C the cross correlation histogram
%         B (optional) a vector with the times corresponding to the bins in msec
%
% example:
%  [C, B] = CrossCorr(t1,t2,2,500);
%  returns a histogram with 2msec bins ranging from -0.5 sec to +0.5 sec centered at zero
%
% version 1.0
