function PI = freqMulti(f)
%FREQMULTI Calculate frequency coefficient for ITD-azimuth warping
%
%   PI = IOSR.AUDITORY.FREQMULTI(F) calculates the coefficient PI for
%   frequency F for use in converting between ITD and azimuth in Kuhn's
%   model [1].
% 
%   References
% 
%   [1] Kuhn, G.F. (1977), Model for the interaural time differences in the
%       azimuthal plane, The Journal of the Acoustical Society of America
%       62, 1, 157-167.
% 
%   See also IOSR.AUDITORY.AZIMUTH2ITD, IOSR.AUDITORY.ITD2AZIMUTH.

%   Copyright 2016 University of Surrey.

    PI = zeros(size(f));
    PI(f<=500) = 3;
    PI(f>=3000) = 2;
    IX = f>500 & f<3000;
    PI(IX) = 2.5+0.5.*cos(pi.*((log2((sqrt(6).*f(IX))./1250))./(log2(6))));

end
