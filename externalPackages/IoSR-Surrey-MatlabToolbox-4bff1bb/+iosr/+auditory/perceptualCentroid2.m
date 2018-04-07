function pc = perceptualCentroid2(x,fs,varargin)
%PERCEPTUALCENTROID2 Alternative perceptual spectral centroid
% 
%   PC = IOSR.AUDITORY.PERCEPTUALCENTROID2(X,FS) calculates the spectral
%   centroid of signal X, sampled at FS, with respect to mel-frequency.
%   
%   The algorithm first calculates the long-term average spectrum (LTAS) of
%   X according to IOSR.DSP.LTAS via the short-time Fourier transform
%   (STFT). The perceptual centroid is then calculated using the LTAS. This
%   is in contrast to IOSR.AUDITORY.PERCEPTUALCENTROID, which calculates
%   the perceptual centroid for a series of spectrogram segments, and then
%   takes the average (and other statistics).
% 
%   X can be a vector, matrix, or multidimensional array;
%   PERCEPTUALCENTROID2 will operate along the first non-signleton
%   dimension, and return a value for each corresponding row/column/etc.
% 
%   PC = IOSR.AUDITORY.PERCEPTUALCENTROID2(X,FS,'PARAMETER',VALUE) allows
%   numerous parameters to be specified. These parameters are:-
%       'cref'     : {27.5} | scalar
%           Specifies the reference frequency when calulating the centroid
%           in terms of cents.
%       'dim'      : {first non-singleton dimension} | integer
%           Specify the dimension over which to calculate the perceptual
%           spectral centroid.
%       'hop' : {512} | integer
%           Specifies the hop size of the STFT.
%       'loudness' : {'none'} | 'A' | 'B' | 'C' | 'D' | 'ISO-226'
%           Specifies whether loudness weighting is applied to the spectrum
%           prior to the centroid calculation. The default is for no
%           weighting to be applied. 'A', 'B', 'C', and 'D'  correspond to
%           frequency weighting curves defined in IEC 61672:2003; 'ISO-226'
%           applies loudness weighting derived from ISO 226:2003.
%       'noct'    : {12} | scalar
%           Apply 1/noct-octave smoothing to the frequency spectrum.
%           Setting 'noct' to 0 results in no smoothing.
%       'output'   : {'units'} | 'hz'
%           Specifies whether the output data are in Hz ('hz') or in units
%           determined by the 'scale' option ('units') (default).
%       'phon'     : {65} | scalar
%           Specifies the loudness level if using ISO 226:2003-based
%           loudness weighting.
%       'scale'    : {'mel'} | 'linear' | 'erb' | 'cents'
%           Specifies frequency scale use to calculate the centroid. 'mel'
%           uses the mel-frequency scale, 'linear' uses a linear frequency
%           scale (corresponding to the traditional spectral centroid
%           measure), 'erb' uses the equivalent-rectangular-bandwidth-rate
%           scale, and cents uses a scale based on musical cents with
%           respect to A-1 (27.5 Hz).
%       'window'   : {hamming(1024)} | vector
%           Specifies the window applied to each STFT segment.
% 
%   Examples
% 
%     Example 1: Calculate the mel-frequency spectral centroid of random
%       noise:
% 
%       fs = 44100;
%       x = randn(fs,1);
%       PC = iosr.auditory.perceptualCentroid2(x,fs);
% 
%     Example 2: Calculate the ERB-spaced spectral centroid of random
%       noise with the output in Hertz:
% 
%       PC = iosr.auditory.perceptualCentroid2(x,fs,...
%           'scale','erb','output','hz');
% 
%     Example 3: Calculate the musically-spaced spectral centroid of
%       random noise using a 4096-point Hann window:
% 
%       PC = iosr.auditory.perceptualCentroid2(x,fs,...
%           'scale','cents','window',hann(4096));
% 
%   Authors: Chris Hummersone & Kirsten Hermes, 2014
%   
%   See also IOSR.DSP.LTAS, IOSR.AUDITORY.PERCEPTUALCENTROID.

%   Copyright 2016 University of Surrey.

    %% read inputs and make default assignments

    assert(~isscalar(x), 'iosr:perceptualCentroid2:invalidX', '''X'' cannot be a scalar')

    dims = size(x);

    propNames = {'dim','hop','noct','window','scale','output','cref','loudness','phon'}; % permitted prop names

    options = struct;
    % read parameter/value options
    if nargin>2 % if parameters are specified
        % count arguments
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('iosr:perceptualCentroid2:nameValuePairs','PERCEPTUALCENTROID2 needs propertyName/propertyValue pairs following the x and fs arguments.')
        end
        % write parameters to options struct
        for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
            if any(strcmpi(pair{1},propNames))
                propName = char(propNames(find(strcmpi(pair{1},propNames),1,'last')));
                options.(propName) = pair{2};
            else
                error('iosr:perceptualCentroid2:unknownOption',['Unknown option ''' pair{1} '''']);
            end
        end
    end

    % other parameters
    dim = getProperty(options,'dim',find(dims>1,1,'first'));
    assert(dim>0 && dim<=length(dims) && isscalar(dim) && round(dim)==dim,...
        'iosr:perceptualCentroid2:invalidDim', ...
        '''dim'' must be greater than zero and less than or equal to the number of dimensions in X.');

    % other parameters
    window = getProperty(options,'window',hamming(1024));
    hop = getProperty(options,'hop',512);
    noct = getProperty(options,'noct',12);
    scale = getProperty(options,'scale','mel');
    output = getProperty(options,'output','units');
    cref = getProperty(options,'cref',27.5);
    loudness = getProperty(options,'loudness','none');

    % tests
    assert(isscalar(cref) && cref>0, 'iosr:perceptualCentroid2:invalidCref', '''cref'' must be greater than 0.');

    %% permute and rehape x to operate down columns

    order = mod(dim-1:dim+length(dims)-2,length(dims))+1;
    dims_shift = dims(order);
    x2 = rearrange(x,order,[dims_shift(1) numel(x)/dims_shift(1)]);

    %% determine functions

    fHandleDoNothing = @(x,~) (x);

    % choose frequency transformation
    switch lower(scale)
        case 'linear'
            fHandleFromF = fHandleDoNothing;
            fHandleToF = fHandleDoNothing;
        case 'mel'
            fHandleFromF = @frequency2mel;
            fHandleToF = @mel2frequency;
        case 'erb'
            fHandleFromF = @frequency2erb;
            fHandleToF = @erb2frequency;
        case 'cents'
            fHandleFromF = @frequency2cents;
            fHandleToF = @cents2frequency;
        otherwise
            error('iosr:perceptualCentroid2:unknownScale',['Requested scale ''' scale ''' not recognised. Options are ''linear'', ''mel'', ''erb'', or ''cents''.']);
    end

    % warn if cref but scale~=cents
    if ~strcmpi(scale,'cents') && ...
            any(cell2mat(strfind(lower(varargin(cellfun(@ischar,varargin))),'cref')))
        warning('iosr:perceptualCentroid2:cref','Option ''cref'' only affects the output when ''scale'' is set to ''cents''.')
    end

    %% calculate centroid and stats
    
    % pre-allocate outputes
    pc = zeros(1,size(x2,2));

    % determine output conversion
    switch lower(output)
        case 'hz'
            % use fHandleToF
        case 'units'
            fHandleToF = fHandleDoNothing;
        otherwise
            error('iosr:perceptualCentroid2:unknownUnits',['Requested output ''' output ''' not recognised. Options are ''hz'', or ''units''.']);
    end

    % determine phon, if relevant
    switch lower(loudness)
        case 'iso-226'
            phon = getProperty(options,'phon',65);
        otherwise
            if isfield(options,'phon')
                warning('iosr:perceptualCentroid2:phon','''phon'' option has no effect unless using ''ISO-226'' loudness weighting.')
            end
    end

    % determine loudness weighting
    switch lower(loudness)
        case 'none'
            w_l = @(x) ones(size(x));
        case 'a'
            w_l = @(x) iosr.auditory.loudWeight(x,'a');
        case 'b'
            w_l = @(x) iosr.auditory.loudWeight(x,'b');
        case 'c'
            w_l = @(x) iosr.auditory.loudWeight(x,'c');
        case 'd'
            w_l = @(x) iosr.auditory.loudWeight(x,'d');
        case 'iso-226'
            w_l = @(x) iosr.auditory.loudWeight(x,phon);
        otherwise
            error('iosr:perceptualCentroid2:unknownLoudness',['Requested loudness weighting ''' loudness ''' not recognised. Options are ''none'', or ''A''.']);
    end

    for c = 1:size(x2,2) % across the dim in input
        % caluclate spectrogram
        [X,f] = iosr.dsp.ltas(x2(:,c),fs,'win',window,'hop',hop,'noct',noct,'units','none');
        % ignore frequencies greater than Nyquist
        IX = f<=fs/2;
        f = f(IX);
        % convert frequencies
        g = fHandleFromF(f,cref);
        % calculate magnitude
        mag = sqrt(X(IX,:));
        % calculate loudness weighting
        loud = w_l(f);
        % calculate centroid
        pc(c) = fHandleToF(sum((loud.*mag).*g)./sum(loud.*mag),cref);
    end

    %% inversely permute output back to input dimensions
    
    pc = irearrange(pc,order,[1 dims_shift(2:end)]);

end

function f = mel2frequency(m,~)
%MEL2FREQUENCY convert mel to frequency
    f = 700.*(exp(m./1127)-1);
end

function m = frequency2mel(f,~)
%FREQUENCY2MEL convert frequency to mel
    m = 1127.*log(1+(f./700));
end

function f = erb2frequency(b,~)
%ERB2FREQUENCY convert ERB to frequency
    f = (10.^(b./21.4)-1)./0.00437;
end

function b = frequency2erb(f,~)
%FREQUENCY2ERB convert frequency to ERB
    b = 21.4.*log10(1+(0.00437.*f));
end

function f = cents2frequency(c,cref)
%CENTS2FREQUENCY convert cents to frequency
    f = cref.*(2.^(c./1200)-eps);
end

function c = frequency2cents(f,cref)
%FREQUENCY2CENTS convert frequency to cents
    c = 1200.*log2((f./cref)+eps);
end

function [propValue,isset] = getProperty(options,propName,default)
%GETPROPERTY return property/default value
    if isfield(options,propName)
        propValue = options.(propName);
        isset = true;
    else
        propValue = default;
        isset = false;
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
