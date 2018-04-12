function w = lapwin(L,b)
%LAPWIN Laplace window.
%   
%   IOSR.DSP.LAPWIN(N) returns an N-point Laplace window. The window w(x) is
%   calculated for -10<=x<=10.
% 
%   IOSR.DSP.LAPWIN(N,B) returns an N-point Laplace window using scale parameter B.
%   If omitted, B is 2.
%   
%   See also GAUSSWIN, CHEBWIN, KAISER, TUKEYWIN, WINDOW.

%   Copyright 2016 University of Surrey.

    assert(isscalar(L) && round(L)==L, 'iosr:lapwin:invalidL', 'L must be a scalar and whole number.')
    
    if nargin<2
        b = 2;
    else
        assert(isscalar(b), 'iosr:lapwin:invalidB', 'b must be a scalar.')
    end
    
    w = zeros(L,1);
    N = L-1;
    x = 10.*((0:N)'-N/2)./(N/2);
    w(x<0) = exp(-(-x(x<0)./b));
    w(x>=0) = exp(-(x(x>=0)./b));

end
