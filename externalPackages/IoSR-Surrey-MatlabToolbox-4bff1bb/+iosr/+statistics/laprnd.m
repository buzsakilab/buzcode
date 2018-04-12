function r = laprnd(varargin)
%LAPRND Pseudorandom numbers drawn from the Laplace distribution
%   
%   R = IOSR.STATISTICS.LAPRND(N) returns an N-by-N matrix containing
%   pseudorandom values drawn from the Laplace distribution with mean = 0
%   and standard deviation = 1.
% 
%   IOSR.STATISTICS.LAPRND(M,N) or IOSR.STATISTICS.LAPRND([M,N]) returns an
%   M-by-N matrix.
%   
%   IOSR.STATISTICS.LAPRND(M,N,P,...) or
%   IOSR.STATISTICS.LAPRND([M,N,P,...]) returns an M-by-N-by-P-by-...
%   array.
% 
%   IOSR.STATISTICS.LAPRND returns a scalar.
% 
%   IOSR.STATISTICS.LAPRND(SIZE(A)) returns an array the same size as A.
%
%   Note: The size inputs M, N, P, ... should be nonnegative integers.
%   Negative integers are treated as 0.
%
%   The sequence of numbers produced by LAPRND is determined by the
%   settings of the uniform random number generator that underlies RAND,
%   RANDI, and RANDN. Control that shared random number generator using
%   RNG.
% 
%   See also IOSR.STATISTICS.TRIRND, RAND, RANDN, RANDI, RNG.

%   Based on code (Matlab FE File ID: #13705) written by Elvis Chen, 2007.

%   Copyright 2016 University of Surrey.

    % Generate Laplacian noise
    u = rand(varargin{:})-0.5;
    b = 1/sqrt(2);
    r = -b*sign(u).*log(1-(2*abs(u)));

end
