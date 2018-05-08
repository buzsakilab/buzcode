function ild = calcIld(L,R,method)
%CALCILD Calculate normalised interaural level difference
% 
%   ILD = IOSR.AUDITORY.CALCILD(L,R) calculates the ILD of vectors L and R.
%   Rather than a logarithmic ratio, the returned ILD is in the range
%   [-1,1]. The algorithm has the following steps:
%   
%   1. Calculate total power for each L and R vectors.
%   2. Calculate difference of powers and divide by power sum.
%   3. Square result, ensuring sign is retained. This improves contrast 
%      when the signals are of a similar level.
% 
%   ILD = IOSR.AUDITORY.CALCILD(L,R,METHOD) determines the nature of the
%   ILD that is returned, and consequently how it is calculated. When
%   method is 'overall' (default), the ILD for the entire signal is
%   returned, based on the power carried by the fine structure. When method
%   is 'vector', the ILD is calculated using the instaneous Hilbert
%   envelope and the function returns a vector the same size as L or R.
% 
%   The function is derived from Ben Supper's thesis "An onset-guided
%   spatial analyser for binaural audio"

%   Copyright 2016 University of Surrey.

    assert(isvector(L) & isvector(R), 'iosr:calcIld:invalidInput', 'L and R must be vectors')
    assert(all(size(L)==size(R)), 'iosr:calcIld:invalidInput', 'L and R must be the same size')

    if nargin<3
        method='overall';
    else
        assert(ischar(method), 'iosr:calcIld:invalidMethod', 'METHOD must be a string.')
    end

    switch lower(method)
        case 'vector'
            envL = calc_env(L);
            envR = calc_env(R);
            p_L = envL.^2;
            p_R = envR.^2;
        case 'overall'
            p_L = sum(L.^2);
            p_R = sum(R.^2);
        otherwise
            error('iosr:calcIld:unknownMethod',['Unknown method ''' method '''.'])
    end
    ild = (p_L-p_R)./(p_L+p_R);
    ild = sign(ild).*(ild.^2);

    end

function y = calc_env(x)
%CALC_ENV return signal Hilbert envelope

    y = abs(hilbert(x));

end
