function rms = rms(x,dim)
%RMS Calculate the rms of a vector or matrix
% 
%   RMS = IOSR.DSP.RMS(X) calculates the Root Mean Square of X along the
%   first non-singleton dimension of vector or matrix X. For vectors,
%   IOSR.DSP.RMS(X) is the RMS value of the elements in x. For matrices,
%   IOSR.DSP.RMS(X) is a row vector containing the RMS value of each
%   column.  For N-D arrays, IOSR.DSP.RMS(X) is the RMS value of the
%   elements along the first non-singleton dimension of X.
% 
%   RMS = IOSR.DSP.RMS(X,DIM) calculates the RMS of X along the dimension
%   DIM.

%   Copyright 2016 University of Surrey.

    if nargin == 1
        dim = find(size(x)~=1,1,'first');
    end

    rms = sqrt(mean(x.^2,dim));

end
