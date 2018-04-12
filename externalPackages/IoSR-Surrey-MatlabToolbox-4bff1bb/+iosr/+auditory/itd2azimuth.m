function azimuth = itd2azimuth(itd,f)
%ITD2AZIMUTH Convert ITD to azimuth
%
%   AZIMUTH = IOSR.AUDITORY.ITD2AZIMUTH(ITD,F) converts the interaural time
%   difference ITD (in seconds) at frequency F (in Hz) to AZIMUTH (in
%   degrees) according to Kuhn's model [1].
% 
%   References
% 
%   [1] Kuhn, G.F. (1977), Model for the interaural time differences in the
%       azimuthal plane, The Journal of the Acoustical Society of America
%       62, 1, 157-167.
% 
%   See also IOSR.AUDITORY.AZIMUTH2ITD, IOSR.AUDITORY.FREQMULTI.

%   Copyright 2016 University of Surrey.

    assert(isnumeric(itd), 'iosr:itd2azimuth:invalidItd', 'ITD must be numeric')
    assert((isscalar(f) | isscalar(itd)) | numel(f)==numel(itd),...
        'iosr:itd2azimuth:invalidInputs', ...
        'F or ITD must be a scalar, or F and ITD must be the same size')

    % Check sanity of ITDs
    assert(all(itd<=iosr.auditory.azimuth2itd(90,f)),...
        'iosr:itd2azimuth:invalidItd', ...
        'ITDs greater than maximum ITD [=iosr.auditory.azimuth2itd(90,f)] have been specified.')

    czero = 344;
    azimuth = asind((itd.*czero)./(iosr.auditory.freqMulti(f).*0.091));

end
