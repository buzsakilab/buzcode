function w = rcoswin(N)
%RCOSWIN Raised cosine window.
%   
%   IOSR.DSP.RCOSWIN(N) returns an N-point raised cosine window.
%   
%   See also GAUSSWIN, CHEBWIN, KAISER, TUKEYWIN, WINDOW.

%   Copyright 2016 University of Surrey.

    assert(isscalar(N) && round(N)==N, 'iosr:rcoswin:invalidN', 'N must be a scalar and whole number.')
    assert(N > 0, 'iosr:rcoswin:invalidN', 'N must be greater than 0.')
    
    switch N
        case 1
            w = 1;
        case 2
            w = [.5; .5];
        otherwise
            r = linspace(-pi,pi,N)';
            w = (cos(r)+1)/2;
    end

end
