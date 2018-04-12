function y = sincFilter(x,Wn,N,dim)
%SINCFILTER Apply a near-ideal low-pass or band-pass brickwall filter
% 
%   Y = IOSR.DSP.SINCFILTER(X,WN) applies a near-ideal low-pass or
%   band-pass brickwall filter to the array X, operating along the first
%   non-singleton dimension (e.g. down the columns of a matrix). The
%   cutoff frequency/frequencies are specified in WN. If WN is a scalar,
%   then WN specifies the low-pass cutoff frequency. If WN is a two-element
%   vector, then WN specifies the band-pass interval. WN must be 0.0 < WN <
%   1.0, with 1.0 corresponding to half the sample rate.
% 
%   The filtering is performed by FFT-based convolution of X with the sinc
%   kernel.
% 
%   Y = IOSR.DSP.SINCFILTER(X,WN,N) allows the filter length to be
%   specified. The default value is N=1025. The filter length is doubled in
%   the band-pass case. In either case, if N is even the final filter
%   length will be N+1.
%   
%   Y = IOSR.DSP.SINCFILTER(X,WN,N,DIM) applies the specified filter along
%   the dimension DIM.
%   
%   Y = IOSR.DSP.SINCFILTER(X,WN,[],DIM) applies the specified filter along
%   the dimension dim using the default filter length.
% 
%   See also IOSR.DSP.CONVFFT.

%   Copyright 2016 University of Surrey.

    %% test input

    assert(nargin>=2, 'iosr:sincFilter:nargin', 'Not enough input arguments')
    assert((numel(Wn)==1 || numel(Wn)==2) & isnumeric(Wn), 'iosr:sincFilter:invalidWn', 'Wn must be a scalar or two-element vector.')
    assert(isnumeric(x), 'iosr:sincFilter:invalidX', 'x must be a numeric array.')
    assert(all(Wn<=1) & all(Wn>=0), 'iosr:sincFilter:invalidWn', 'Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to half the sample rate.')
    dims = size(x);
    if nargin<4
        dim = [];
    else
        assert(isint(dim) & numel(dim)==1, 'iosr:sincFilter:invalidDim', 'dim must be an integer')
        assert(dim<=length(dims), 'iosr:sincFilter:invalidDim', 'dim must be less than or equal to the number of dimensions in x')
    end
    if nargin<3
        N = [];
    elseif ~isempty(N)
        assert(isscalar(N) & isint(N), 'iosr:sincFilter:invalidN', 'N must be an integer scalar.')
    end

    %% assign defaults

    if isempty(N)
        N = 1025;
    end
    if isempty(dim)
        dim = find(dims>1,1,'first');
    end
%%


%take gpuArray as input
d = whos('x');
useGPU = false;
if strcmp(d.class,'gpuArray')
    useGPU = true;
end



    %% reshape input to matrix with requested dim made as first dimension

    % reshape data so function works down columns
    order = mod(dim-1:dim+length(dims)-2,length(dims))+1;
    x = permute(x,order);
    dims_shift = dims(order);
    x = reshape(x,dims_shift(1),numel(x)/dims_shift(1));
    y = zeros(size(x));
    

    %% create filter kernel

    if numel(Wn)==1 % low-pass
        n = -floor(N/2):floor(N/2); % create time base
        B = sinc_kernel(Wn,n); % make kernel
    else % band-pass
        n = -N:N; % create time base
        B = sinc_kernel(Wn(2),n)-sinc_kernel(Wn(1),n); % make kernel
    end

    %% apply filter

    
    %%
    if useGPU
        y = gpuArray(y);
        B = gpuArray(B);
        
    end
    
    
    %%
    for N = 1:size(y,2)
        y(:,N) = iosr.dsp.convFft(x(:,N),B,'same');
    end

    %% reshape out to match input

    y = reshape(y,dims_shift);
    y = ipermute(y,order);

end

function k = sinc_kernel(Wn,n)
% SINC_KERNEL Make sinc kernel

    k = sinc(Wn*n).*Wn;

end

function y = isint(x)
% ISINT Test whether x is integer value (not type)

    y = x==round(x);

end
