function y=erbRate2hz(x)
%ERBRATE2HZ Convert ERB rate to Hz.
% 
%   Y = IOSR.AUDITORY.ERBRATE2HZ(X) converts the ERB number X to the
%   eqivalent frequency Y (in Hz).
% 
% See also IOSR.AUDITORY.HZ2ERBRATE.

%   Copyright 2016 University of Surrey.

    y=(10.^(x/21.4)-1)/4.37e-3;

end
