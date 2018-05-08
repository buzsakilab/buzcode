function [s,f] = ltas(x,fs,varargin)
%LTAS calculate the long-term average spectrum of a signal
% 
%   S = IOSR.DSP.LTAS(X,FS) calculates the long-term average spectrum
%   (LTAS) of signal X, sampled at FS Hz. The spectrum is calculated from
%   the average power spectral density (PSD) obtained from a series of
%   overlapping FFTs; the FFT length is 4096, and the hop size is 2048. The
%   segments of X are Hann-windowed. The average PSD is then
%   Gaussian-smoothed to 1/3-octave resolution.
% 
%   X can be a vector, matrix, or multidimensional array; LTAS will operate
%   along the first non-signleton dimension, and return the LTAS for each
%   corresponding row/column/etc.
%   
%   S = IOSR.DSP.LTAS(X,FS,'PARAMETER','VALUE') allows numerous parameters
%   to be specified. These parameters are:-
%       'dim'     : {find(size(X)>1,1,'first')} | scalar
%           Specifies the dimension of operation (defaults to the first
%           non-singleton dimension).
%       'graph'   : {false} | true
%           Choose whether to plot a graph of the LTAS.
%       'hop'     : {NFFT/2} | scalar
%           Specifies the step size through X used to calculate each
%           segment. NFFT is determined by the 'win' parameter.
%       'noct'    : {3} | scalar
%           Apply 1/noct-octave smoothing to the frequency spectrum.
%           Setting 'noct' to 0 results in no smoothing.
%       'scaling' : {'none'} | 'max0'
%           Specifies any scaling to apply to S. By default, no scaling is
%           applied. If scaling is set to 'max0', S will be scaled to have
%           a maximum value of 0dB.
%       'units'   : {dB} | 'none'
%           Specifies the output units. By default the PSD is calculated in
%           dB. Otherwise the PSD is returned directly.
%       'win'     : {4096} | scalar | vector
%           Specifies the window or FFT length NFFT used to calculate the
%           spectrum. If 'win' is a scalar, it specifies the FFT length,
%           and a Hann window is applied to each segment. If 'win' is a
%           vector, NFFT is the length of the vector, and the vector is
%           multiplied with each segment.
% 
%   [S,F] = IOSR.DSP.LTAS(...) returns the frequencies F for each
%   corresponding bin of S.
% 
%   Example
%   
%       % Plot the 1/6th-octave-smoothed LTAS of the Handel example
%       load handel.mat
%       figure
%       iosr.dsp.ltas(y,Fs,'noct',6,'graph',true);
%   
%   See also IOSR.DSP.STFT, IOSR.DSP.SMOOTHSPECTRUM.

%   Copyright 2016 University of Surrey.

    %% parse input
    
    if nargin < 2
        error('iosr:ltas:nargin''Not enough input arguments')
    end
    
    options = struct(...
        'win',4096,...
        'hop',[],...
        'noct',3,...
        'dim',[],...
        'units','db',...
        'scaling','none',...
        'graph',false);

    % read parameter/value inputs
    if nargin>2 % if parameters are specified
        % read the acceptable names
        optionNames = fieldnames(options);
        % count arguments
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('iosr:ltas:nameValuePair','LTAS needs propertyName/propertyValue pairs')
        end
        % overwrite defults
        for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
           IX = strcmpi(pair{1},optionNames); % find match parameter names
           if any(IX)
              % do the overwrite
              options.(optionNames{IX}) = pair{2};
           else
              error('iosr:ltas:unknownOption','%s is not a recognized parameter name',pair{1})
           end
        end
    end
    
    %% check and assign input
    
    % required inputs
    assert(isnumeric(x), 'iosr:ltas:invalidX', 'X must be numeric.');
    assert(isscalar(fs), 'iosr:ltas:invalidFs', 'FS must be a scalar.');
    assert(fs>0, 'iosr:ltas:invalidFs', 'FS must be greater than 0.');
    assert(isint(fs), 'iosr:ltas:invalidFs', 'FS must be an integer.');
    
    % determine fft parameters
    if numel(options.win)>1 && isvector(options.win)
        NFFT = length(options.win);
        win = options.win;
    elseif isscalar(options.win)
        NFFT = options.win;
        win = NFFT;
    else
        error('iosr:ltas:invalidWin','''WIN'' must be a vector or a scalar');
    end
    assert(isint(NFFT) && NFFT>0, 'iosr:ltas:invalidNfft', '''WIN'' must be a positive integer');
    
    % determine hop
    hop = options.hop;
    if isempty(hop)
        hop = fix(NFFT/2);
    end
    
    % smoothing
    Noct = options.noct;
    
    % dimension to operate along
    dim = options.dim;
    dims = size(x);
    if isempty(dim)
        dim = find(dims>1,1,'first');
    else
        assert(isnumeric(dim),'iosr:ltas:invalidDim', 'DIM must be an integer');
        assert(isint(dim), 'iosr:ltas:invalidDim', 'DIM must be an integer or empty');
        assert(dim>0, 'iosr:ltas:invalidDim', 'DIM must be greater than 0')
    end
    
    %% permute and rehape x to operate down columns
    
    % number of useful coefficients
    if mod(NFFT,2)==0
        Nout = (NFFT/2)+1;
    else
        Nout = (NFFT+1)/2;
    end
    
    % reorder and permute
    order = mod(dim-1:dim+length(dims)-2,length(dims))+1;
    dims_shift = dims(order);
    x = rearrange(x,order,[dims_shift(1) prod(dims_shift(2:end))]);
    dims_out_shift = dims_shift;
    dims_out_shift(1) = Nout;
    
    %% calculate spectra
    
    % choose units function
    switch lower(options.units)
        case 'db'
            units = @(x) 10*log10(x);
            max0scale = @(s) s-max(s(:));
            labelY = 'Power spectral density [dBFS]';
            yscale = 'linear';
        case 'none'
            units = @(x) x;
            max0scale = @(s) s./max(s(:));
            labelY = 'Power spectral density';
            yscale = 'log';
        otherwise
            error('iosr:ltas:unknownUnits','Unknown units option ''%s''',options.units);
    end
    
    % do calculations
    s = zeros(dims_out_shift);
    for c = 1:size(x,2)
        [S,f] = iosr.dsp.stft(x(:,c),win,hop,fs); % short-time ft
        s(:,c) = mean(abs(S/NFFT).^2,2); % mean PSD
        s(:,c) = units(s(:,c)); % put into dB
        s(:,c) = iosr.dsp.smoothSpectrum(s(:,c),f,Noct); % smooth
    end
    
    % invert permutation
    s = irearrange(s,order,dims_out_shift);
    
    % scale output
    switch lower(options.scaling)
        case 'max0'
            s = max0scale(s);
        case 'none'
            % do nothing
        otherwise
            error('iosr:matchEQ:unknownScaling','Unknown scaling option ''%s''',options.scaling);
    end
    
    %% plot
    
    assert(islogical(options.graph) && numel(options.graph)==1, 'iosr:ltas:invalidGraph', '''graph'' option must be logical.')
    if options.graph
        semilogx(f,rearrange(s,order,[dims_out_shift(1) prod(dims_out_shift(2:end))]));
        xlabel('Frequency [Hz]');
        ylabel(labelY);
        set(gca,'yscale',yscale);
        grid on
        if prod(dims_out_shift(2:end)) > 1
            legend(num2str((1:prod(dims_out_shift(2:end)))'));
        end
    end
    
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
