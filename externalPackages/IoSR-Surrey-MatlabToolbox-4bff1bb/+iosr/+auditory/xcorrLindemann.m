function [c,lags] = xcorrLindemann(L,R,fs,maxlag,dim)
%XCORRLINDEMANN Cross-correlation based on Lindemann's precedence model
% 
%   C = IOSR.AUDITORY.XCORRLINDEMANN(L,R,FS) calculates the
%   cross-correlation of vectors L and R, sampled at FS Hz, using a maximum
%   cross-correlation lag of 0.001*fs samples (1 ms). The cross-correlation
%   function is based on Lindemann's [1,2] precedence model. C is the same
%   size as L or R, except that size(C,1) = 2*0.001*fs+1.
% 
%   C = IOSR.AUDITORY.XCORRLINDEMANN(L,R,FS,MAXLAG) calculates the
%   cross-correlation using a maximum cross-correlation lag of MAXLAG
%   samples.
% 
%   C = IOSR.AUDITORY.XCORRLINDEMANN(L,R,FS,MAXLAG,DIM) performs the
%   cross-correlation along the dimension DIM. C is the same size as L or
%   R, except that size(C,DIM) = 2*MAXLAG+1.
% 
%   [C,LAGS] = IOSR.AUDITORY.XCORRLINDEMANN(...) returns the lags, in
%   seconds, over which the cross-correlation was calculated.
% 
%   References
%   
%   [1] Lindemann, W. (1986), Extension of a binaural cross-correlation
%     model by contralateral inhibition. I. Simulation of lateralization
%     for stationary signals, The Journal of the Acoustical Society of
%     America 80, 6, 1608-1622.
% 
%   [2] Lindemann, W. (1986), Extension of a binaural cross-correlation
%     model by contralateral inhibition. II. The law of the first wave
%     front, The Journal of the Acoustical Society of America 80, 6,
%     1623-1630.
%   
%   See also IOSR.AUDITORY.LINDEMANNINH, XCORR.

%   Copyright 2016 University of Surrey.

    %% check input
    assert(isequal(size(L),size(R)), 'iosr:xcorrLindemann:invalidSignals', 'L and R must be the same size');
    assert(isscalar(fs) & isnumeric(fs), 'iosr:xcorrLindemann:invalidFs', 'FS must be a scalar');
    
    % check for maxlag
    if nargin<4
        maxlag = 0.001*fs;
    else
        assert(isscalar(maxlag) & isnumeric(maxlag), 'iosr:xcorrLindemann:invalidMaxlag', 'MAXLAG must be a scalar');
    end
    
    % check for dim
    dims = size(L);
    if nargin<5
        dim = find(dims>1,1,'first');
    else
        assert(isscalar(dim) && round(dim)==dim, 'iosr:xcorrLindemann:invalidDim', 'DIM must be an integer scalar')
    end

    % check L and R have a valid range
    if any(L(:)<0) || any(R(:)<0) || any(L(:)>1) || any(R(:)>1)
        error('iosr:xcorrLindemann:inputOutOfRange','L and/or R contain values outside of the range [0,1]. This is not allowed. Use LINDEMANN_INH to pre-process the L and R inputs.')
    end

    % check C function is compiled
    iosr.general.checkMexCompiled('-largeArrayDims',fullfile(fileparts(mfilename('fullpath')),'xcorrLindemann_c.c'))
    
    %% re-arrange input

    % re-arrange array to operate along columns
    order = mod(dim-1:dim+length(dims)-2,length(dims))+1;
    dims_shift = dims(order);
    L = rearrange(L,order,[dims_shift(1),numel(L)/dims_shift(1)]);
    R = rearrange(R,order,[dims_shift(1),numel(R)/dims_shift(1)]);
    
    %% do cross-correlation
    
    % time constants
    t_int = 5; % ms
    t_inh = 10; % ms
    
    xcorr_length = round(2*maxlag)+1; % length of cross-correlation
    c = zeros(xcorr_length,size(L,2)); % pre-allocate
    for n = 1:size(L,2)
        c(:,n) = iosr.auditory.xcorrLindemann_c(L(:,n),R(:,n),fs,maxlag,t_inh,t_int);
    end
    
    %% return values
    
    % rearrange to match input dimensions
    c = irearrange(c,order,[xcorr_length dims_shift(2:end)]);

    % return lags
    if nargin>1
        lags = (-maxlag:maxlag)./fs;
    end

end

function y = rearrange(x,dim_order,shape)
%REARRANGE rearrange data to 2-D matrix with target dim as column
    
    y = permute(x,dim_order);
    y = reshape(y,shape);
    
end

function y = irearrange(x,dim_order,shape)
%IREARRANGE inverse rearrangement
    
    y = reshape(x,shape);
    y = ipermute(y,dim_order);
    
end
