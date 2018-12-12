function [data] = bz_shuffleCellID(data)
% randomly permutes cell order

data = data(randperm(size(data,1)),:);