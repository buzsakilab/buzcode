function itd = azimuth2itd(azimuth,f)
%AZIMUTH2ITD Convert azimuth in degrees to ITD
% 
%   ITD = IOSR.AUDITORY.AZIMUTH2ITD(AZIMUTH,F) converts the azimuth AZIMUTH
%   in degrees at frequency F to interaural time difference according to
%   Kuhn's model [1].
% 
%   References
% 
%   [1] Kuhn, G.F. (1977), Model for the interaural time differences in the
%       azimuthal plane, The Journal of the Acoustical Society of America
%       62, 1, 157-167.
% 
%   See also IOSR.AUDITORY.ITD2AZIMUTH, IOSR.AUDITORY.FREQMULTI.

%   Copyright 2016 University of Surrey.

    assert(isnumeric(azimuth), 'iosr:azimuth2itd:invalidInput', 'AZIMUTH must be numeric')
    assert((isscalar(f) | isscalar(azimuth)) | numel(f)==numel(azimuth),...
        'iosr:azimuth2itd:invalidInput', ...
        'F or ITD must be a scalar, or F and AZIMUTH must be the same size')

    czero = 344;
    itd = (iosr.auditory.freqMulti(f).*0.091.*sind(azimuth))./czero;

end
