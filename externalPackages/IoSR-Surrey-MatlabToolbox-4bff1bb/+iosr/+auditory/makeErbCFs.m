function cfs = makeErbCFs(mincf,maxcf,numchans)
%MAKEERBCFS Make a series of center frequencies equally spaced in ERB-rate.
% 
%   This function makes a vector of center frequenies equally spaced on the
%   ERB-rate scale.
%   
%   CFS = IOSR.AUDITORY.MAKEERBCFS(MINCF,MAXCF,NUMCHANS) creates NUMCHANS
%   centre frequencies between MINCF and MAXCF.
%
%   Adapted from code written by: Guy Brown, University of Sheffield, and
%   Martin Cooke.
% 
%   See also IOSR.AUDITORY.ERBRATE2HZ, IOSR.AUDITORY.HZ2ERBRATE.

%   Copyright 2016 University of Surrey.

    cfs = iosr.auditory.erbRate2hz(...
        linspace(iosr.auditory.hz2erbRate(mincf),...
        iosr.auditory.hz2erbRate(maxcf),numchans));

end
