function y=hz2erbRate(x)
%HZ2ERBRATE Convert Hz to ERB rate
%
%   Y = IOSR.AUDITORY.HZ2ERBRATE(X) converts the frequency X (in Hz) to the
%   eqivalent ERB number Y.
% 
%   See also IOSR.AUDITORY.ERBRATE2HZ, IOSR.AUDITORY.MAKEERBCFS.

%   Copyright 2016 University of Surrey.

    y=(21.4*log10(4.37e-3*x+1));

end
