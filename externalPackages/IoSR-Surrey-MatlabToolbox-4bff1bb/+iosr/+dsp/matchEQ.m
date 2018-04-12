function [y,b] = matchEQ(x,fs,mag,f,varargin)
%MATCHEQ Match the LTAS of a signal to an arbitrary spectral magnitude
% 
%   Y = IOSR.DSP.MATCHEQ(X,FS,MAG,F) matches the long-term average spectrum
%   of X, sampled at FS Hz, to the spectral magnitudes contained in vector
%   MAG sampled at frequencies in vector F. X can be a vector, matrix, or
%   multidimensional array; MATCHEQ will operate along the first
%   non-signleton dimension. Y is the same size as X.
% 
%   Y = IOSR.DSP.MATCHEQ(X,FS,MAG,F,'PARAMETER','VALUE') allows numerous
%   parameters to be specified. These parameters are:-
%       'dim'        : {find(size(X)>1,1,'first')} | scalar
%           Specifies the dimension of operation (defaults to the first
%           non-singleton dimension).
%       'boostRatio' : {1} | scalar
%           Determine the correction ratio for spectral regions that are
%           amplified. By default, the 'ratio' parameter is used. See
%           'ratio' below.
%       'cutRatio'   : {1} | scalar
%           Determine the correction ratio for spectral regions that are
%           attenuated. By default, the 'ratio' parameter is used. See
%           'ratio' below.
%       'hop'        : {NFFT/2} | scalar
%           Specifies the step size through X used to calculate each
%           segment for the LTAS. NFFT is determined by the 'win'
%           parameter.
%       'ltasNorm'   : {'none'} | 'peak' | 'sum'
%           Apply normalisation to the LTAS. By default, no normalisation
%           is applied, and the equalisation will attempt to match the
%           targtet magnitude. By setting the normalisation to 'peak', MAG
%           and the LTAS of X are normalised in order to have a peak gain
%           of unity before calculating the correction filter. By setting
%           the normalisation to 'sum', MAG and the LTAS of X are
%           normalised in order to sum to unity before calculating the
%           correction filter.
%       'noct'       : {3} | scalar
%           Apply 1/noct-octave smoothing to the LTAS. Setting 'noct' to 0
%           results in no smoothing.
%       'order'      : {2048} | scalar
%           Determine the order of the correction filter.
%       'ratio'      : {1} | scalar
%           Determine the correction ratio. A value of 0 means that no
%           correction is applied; 1 means that the full correction is
%           applied; a value greater than 1 will over-correct.
%       'win'        : {4096} | scalar | vector
%           Specifies the window or FFT length NFFT used to calculate the
%           LTAS. If 'win' is a scalar, it specifies the FFT length,
%           and a Hann window is applied to each segment. If 'win' is a
%           vector, NFFT is the length of the vector, and the vector is
%           multiplied with each segment.
% 
%   [Y,B] = IOSR.DSP.MATCHEQ(...) returns the filter coefficients to the
%   array B. B is the same size as X, except that the DIMth dimension of B
%   has length N+1, where N is the filter order and DIM is the dimension of
%   operation.
% 
%   Example
%   
%       % flatten the spectrum of the Handel example
%       load handel.mat
%       f = linspace(0,1,1024);
%       mag = ones(size(f));
%       z = iosr.dsp.matchEQ(y,Fs,mag,f,'ltasNorm','sum');
% 
%   See also IOSR.DSP.LTAS.

%   Copyright 2016 University of Surrey.

    %% parse input
    
    if nargin < 4
        error('iosr:matchEQ:nargin','Not enough input arguments')
    end
    
    options = struct(...
        'win',4096,...
        'hop',[],...
        'noct',3,...
        'dim',[],...
        'ratio',1,...
        'boostRatio',[],...
        'cutRatio',[],...
        'order',2048,...
        'ltasNorm','none');

    % read parameter/value inputs
    if nargin>4 % if parameters are specified
        % read the acceptable names
        optionNames = fieldnames(options);
        % count arguments
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('iosr:matchEQ:nameValuePair','MATCHEQ needs propertyName/propertyValue pairs')
        end
        % overwrite defults
        for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
           IX = strcmpi(pair{1},optionNames); % find match parameter names
           if any(IX)
              % do the overwrite
              options.(optionNames{IX}) = pair{2};
           else
              error('iosr:matchEQ:unknownOption','%s is not a recognized parameter name',pair{1})
           end
        end
    end
    
    %% check and assign input
    
    % required inputs
    assert(isnumeric(x), 'iosr:matchEQ:invalidX', 'X must be numeric.');
    assert(isscalar(fs), 'iosr:matchEQ:invalidFs', 'FS must be a scalar.');
    assert(fs>0, 'iosr:matchEQ:invalidFs', 'FS must be greater than 0.');
    assert(isint(fs), 'iosr:matchEQ:invalidFs', 'FS must be an integer.');
    assert(isnumeric(mag) && isvector(mag), 'iosr:matchEQ:invalidMag', 'MAG must be a numeric vector.')
    assert(isnumeric(f) && isvector(f), 'iosr:matchEQ:invalidF', 'F must be a numeric vector.')
    assert(isequal(size(mag),size(f)), 'iosr:matchEQ:invalidInputs', 'F and MAG must be the same size.')
    assert(all(f>=0), 'iosr:matchEQ:invalidF', 'Frequencies in F must be greater than or equal to 0.')
    assert(all(mag>=0), 'iosr:matchEQ:invalidMag', 'Magnitudes in MAG must be greater than or equal to 0.')
    
    % dimension to operate along
    dim = options.dim;
    dims = size(x);
    if isempty(dim)
        dim = find(dims>1,1,'first');
    else
        assert(isnumeric(dim), 'iosr:matchEQ:invalidDim', 'DIM must be an integer');
        assert(isint(dim), 'iosr:matchEQ:invalidDim', 'DIM must be an integer or empty');
        assert(dim>0, 'iosr:matchEQ:invalidDim', 'DIM must be greater than 0')
    end
    
    f = f(:);
    mag = mag(:);
    
    %% Determine ratios
    
    assert(options.ratio>=0 && isscalar(options.ratio), 'iosr:matchEQ:invalidRatio', 'Ratio must be a positive scalar.')
    
    if isempty(options.boostRatio)
        boostRatio = options.ratio;
    else
        boostRatio = options.boostRatio;
    end
    
    if isempty(options.cutRatio)
        cutRatio = options.ratio;
    else
        cutRatio = options.cutRatio;
    end
    
    assert(boostRatio>=0 && isscalar(boostRatio), 'iosr:matchEQ:invalidBoostRatio', 'BoostRatio must be a positive scalar.')
    assert(cutRatio>=0 && isscalar(cutRatio), 'iosr:matchEQ:invalidCutRatio', 'CutRatio must be a positive scalar.')
    
    %% Determine normalisation
    
    switch lower(options.ltasNorm)
        case 'none'
            fltas = @(x) x;
        case 'sum'
            fltas = @(x) x./sum(x);
        case 'peak'
            fltas = @(x) x./max(x);
        otherwise
            error('iosr:matchEQ:unknownNorm','Unknown normalisation mode ''%s''',options.ltasNorm)
    end
    
    %% permute and rehape x to operate down columns
    
    % reorder and permute
    order = mod(dim-1:dim+length(dims)-2,length(dims))+1;
    dims_shift = dims(order);
    x = rearrange(x,order,[dims_shift(1) prod(dims_shift(2:end))]);
    dims_out_shift = dims_shift;
    b_out_dims = dims_shift; b_out_dims(1) = options.order+1;
    
    %% EQ
    
    assert(isscalar(options.order) && options.order>0 && isint(options.order),...
        'iosr:matchEQ:invalidOrder', ...
        'ORDER must be a positive scalar integer.')
    
    % do calculations
    y = zeros(size(x));
    b = zeros(b_out_dims);
    for c = 1:size(x,2)
        [mX,fX] = iosr.dsp.ltas(x(:,c),fs,...
            'hop',options.hop,...
            'noct',options.noct,...
            'win',options.win,...
            'units','none');
        mX = fltas(mX);
        if ~isequal(fX(:),f)
            mag = interp1(f,mag,fX,'spline','extrap');
        end
        mag = fltas(mag);
        matcher = sqrt(mag./mX);
        matcher(matcher>1) = matcher(matcher>1).*boostRatio;
        matcher(matcher<1) = matcher(matcher<1).*cutRatio;
        b(:,c) = fir2(options.order,2.*fX./fs,matcher);
        y(:,c) = filter(b(:,c),1,x(:,c));
    end
    
    % invert permutation
    y = irearrange(y,order,dims_out_shift);
    b = irearrange(b,order,b_out_dims);
    
end

function y = rearrange(x,order,shape)
%REARRANGE reshape and permute to make target dim column
    y = permute(x,order);
    y = reshape(y,shape);
end

function y = irearrange(x,order,shape)
%IREARRANGE reshape and permute to original size
    y = reshape(x,shape);
    y = ipermute(y,order);
end

function y = isint(x)
%ISINT check if input is whole number
    y = x==round(x);
end
