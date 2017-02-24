function [signals, PC, V, data] = pca1(data)
% PCA1: Perform PCA using covariance.
%
%   Usage:
%       [signals, PC, V] = pca1(data)
%
%   data (in)   - MxN matrix of input data
%             (M dimensions, N trials)
%   signals - MxN matrix of projected data
%   PC      - each column is a principal component
%   V       - Mx1 matrix of variances
%   data (out) - mean subtracted data.

%
% From Jon Schlens PCA tutorial.

[M, N] = size(data);

% subtract off the mean for each dimension
mn = mean(data,2);
data = data- repmat(mn,1,N);

% calculate the covariance matrix
covariance = 1/(N-1)*data*data';

if (isnan(covariance))
    covariance(:,:) = 0
end

% find the eigenvectors and eigenvalues
[PC, V] = eig(covariance);

% extract diagonal of matrix as vector
V = diag(V);

% sort the variances in decreasing order
[junk, rindices] = sort(-1*V);
V = V(rindices);
PC = PC(:, rindices);

% project the original data set
signals = PC'*data;