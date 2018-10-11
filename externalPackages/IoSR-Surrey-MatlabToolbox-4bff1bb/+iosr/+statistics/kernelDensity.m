function [d, xd, bw] = kernelDensity(x, bins, h, kernel)
%KERNELDENSITY Calculate the kernel density of a data set
%   
%   [D, XD] = IOSR.STATISTICS.KERNELDENSITY(X) calculate the kernel density
%   D of a dataset X for query points XD. The kernel density is calculated
%   for 100 query points equally spaced between the minimum and maximum of
%   the data X. The density is estimated using a gaussian kernel with a
%   width that is optimal for normal data. X may be a vector, matrix, or
%   multi-dimensional array; the entire array is treated as the sample. D
%   and XD are 100-point column vectors. NaN are excluded from the
%   calculations.
% 
%   ... = IOSR.STATISTICS.KERNELDENSITY(X, BINS) calculates the density for
%   the query points specified by BINS. If BINS is a scalar, then BINS
%   points are queried between MIN(X(:)) and MAX(X(:)); if BINS is a
%   vector, then the values are used as the query points directly. D and XD
%   are column vectors.
% 
%   ... = IOSR.STATISTICS.KERNELDENSITY(X, BINS, H) uses the bandwidth H to
%   calculate the kernel density. H must be a scalar. BINS may be an empty
%   array in order to use the default described above.
% 
%   ... = IOSR.STATISTICS.KERNELDENSITY(X, BINS, [], KERNEL) uses the
%   kernel function specified by KERNEL to calculate the density. The
%   kernel may be:
%       - 'normal' (default),
%       - 'uniform',
%       - 'triangular',
%       - 'epanechnikov',
%       - 'quartic',
%       - 'triweight',
%       - 'tricube',
%       - 'cosine',
%       - 'logistic',
%       - 'sigmoid', or
%       - 'silverman'.
%   For the uniform case the bandwidth is set to 15% of the range of the
%   data [1]. Otherwise the bandwidth is chosen to be optimal for normal
%   data assuming a gaussian kernel.
% 
%   ... = IOSR.STATISTICS.KERNELDENSITY(X, BINS, H, KERNEL) allows the
%   bins BINS and bandwidth H to be specified directly.
% 
%   [D, XD, BW] = IOSR.STATISTICS.KERNELDENSITY(...) returns the badwidth
%   BW.
% 
%   Examples
% 
%     Example 1: Plot the kernel density of gaussian data
%       figure
%       % gaussian random numbers
%       y = randn(100000, 1);
%       % density
%       [d, xd] = iosr.statistics.kernelDensity(y);
%       % plot
%       plot(xd, d);
% 
%     Example 2: Density trace with 200 bins of width of 10% of data range
%       figure
%       % random numbers
%       y = randn(100000, 1);
%       % y range
%       range = max(y(:)) - min(y(:));
%       % density trace
%       [d, xd] = iosr.statistics.kernelDensity(y, 200, 0.1*range, 'uniform');
%       % plot
%       plot(xd, d);
% 
%   References
% 
%   [1] Hintze, Jerry L.; Nelson, Ray D. (1998). "Violin Plots: A Box
%       Plot-Density Trace Synergism". The American Statistician. 52 (2):
%       181?4. 

    %% input check

    x = x(:);
    x = x(~isnan(x));
    assert(numel(x) > 1, 'X must be a vector, matrix, or array.')
    
    % x bins
    if nargin < 2
        bins = [];
    end
    if isempty(bins)
        bins = 100;
    end
    if isscalar(bins)
        bins = linspace(min(x), max(x), round(bins));
    end
    bins = bins(:);
    
    % bin width
    if nargin < 3
        h = [];
    end
    
    % kernel
    if nargin < 4
        kernel = [];
    end
    if isempty(kernel)
        kernel = 'normal';
    end
    % return kernel function
    switch lower(kernel)
        case 'uniform'
            K = @(u) 0.5*(abs(u) <= 1);
            if isempty(h)
                h = 0.15 * (max(x) - min(x));
            end
        case 'normal'
            K = @(u) ((1/sqrt(2*pi)) * exp(-0.5*(u.^2)));
        case 'triangular'
            K = @(u) ((1-abs(u)) .* (abs(u) <= 1));
        case 'epanechnikov'
            K = @(u) ((0.75*(1-(u.^2))) .* (abs(u) <= 1));
        case 'quartic'
            K = @(u) (((15/16)*(1-(u.^2)).^2) .* (abs(u) <= 1));
        case 'triweight'
            K = @(u) (((35/32)*(1-(u.^2)).^3) .* (abs(u) <= 1));
        case 'tricube'
            K = @(u) (((70/81)*(1-(abs(u).^3)).^3) .* (abs(u) <= 1));
        case 'cosine'
            K = @(u) (((pi/4)*cos((pi/2)*u)) .* (abs(u) <= 1));
        case 'logistic'
            K = @(u) (1 / (exp(u) + 2 + exp(-u)));
        case 'sigmoid'
            K = @(u) ((2/pi) * (1 / (exp(u) + exp(-u))));
        case 'silverman'
            K = @(u) (0.5 * exp((-abs(u))/(sqrt(2))) .* sin((abs(u))/(sqrt(2)) + (pi/4)));
        otherwise
            error('Unknown kernel specified');
    end
    if isempty(h)
        h = ((4*(std(x).^5))/(3*numel(x))).^(1/5);
    end
    
    assert(isscalar(h), 'h must be a scalar')
    
    %% calculate kernel density

    xd = sort(bins);
    d = zeros(size(xd));
    
    for i = 1:numel(xd)
        d(i) = sum(K((x-xd(i))/h))./(numel(x)*h);
    end
    d(isnan(d)) = 0;
    bw = h;

end