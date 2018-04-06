function output = resynthesise(x,fs,cfs,m,varargin)
%RESYNTHESISE Resynthesise a target from a time-frequency mask
% 
%   OUTPUT = IOSR.BSS.RESYNTHESISE(X,FS,CFS,M) takes an input vector X,
%   sampled at FS Hz, and applies the time-frequency mask M. The mask M is
%   a matrix, with one column for each frequency channel, and one row for
%   each sample. A row may also correspond to a frame - see 'frame_length'
%   option below. The centre frequency of each column is contained in CFS.
% 
%   OUTPUT = IOSR.BSS.RESYNTHESISE(X,FS,CFS,M,'PARAMETER',VALUE) allows a
%   number of options to be specified. The options are:
% 
%   %   ({} indicates the default value)
% 
%   'frame_length' : {1} | scalar
%       The frame length for each time-frequency unit in M, in samples.
%   'delay'        : {zeros(size(m,2),1)} | vector
%       The delay in the filterbank when the time-frequency mask M was
%       calculated. The delay is removed only when 'frame_length' > 1 and
%       the full mask must be reconstructed. The value should be a vector
%       whose length is equal to the number of columns in M; delay(n) is
%       the delay in frequency band cfs(n). The delay does not affect the
%       ouput of the reconstruction filters, since the filters are
%       zero-phase.
%   'kernel'       : {1} | vector
%       Allows smoothing to be applied to the full mask. See
%       IOSR.BSS.GETFULLMASK for more details.
%   'filter'       : {'gammatone'} | 'sinc'
%       Allows the reconstruction filter to be specified. By default a
%       gammatone filterbank is used, and cfs corresponds to the filter
%       centre frequencies. Alternatively a series sinc band-pass filters
%       can be used. In this case it is assumed that the frequencies in cfs
%       correspond to the upper cut-off frequencies of the sinc filters.
%       The function IOSR.BSS.CFS2FCS can be used to calculate the cut-off
%       frequencies.
%   'order'        : {fs} | scalar
%       The sinc filter order.
% 
%   Algorithm
% 
%   The algorithm has the following steps:
%   1) If necessary, obtain the full (smoothed) binary mask. See
%      IOSR.BSS.GETFULLMASK.
%   2) Pass the vector x through the reconstruction filterbank, with
%      centre/cut-off frequencies cfs.
%   3) Multiply the full mask with the output of the filterbank.
%   4) Sum the responses of the filterbank to produce a vector the same
%      size as x.
% 
%   See also IOSR.BSS.CFS2FCS, IOSR.BSS.GETFULLMASK, IOSR.BSS.SINCFILTER.

%   Copyright 2016 University of Surrey.

    %% parse input

    if nargin<4
        error('iosr:resynthesise:tooFewInputArgs','Not enough input arguments')
    end

    numchans = size(m,2);

    options = struct(...
        'frame_length',1,...
        'delay',zeros(numchans,1),...
        'kernel',1,...
        'filter','gammatone',...
        'order',fs);

    % read parameter/value inputs
    if nargin > 4 % if parameters are specified
        % read the acceptable names
        optionNames = fieldnames(options);
        % count arguments
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('iosr:resynthesise:nameValuePairs','RESYNTHESISE needs propertyName/propertyValue pairs')
        end
        % overwrite defults
        for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
           IX = strcmpi(pair{1},optionNames); % find match parameter names
           if any(IX)
              % do the overwrite
              options.(optionNames{IX}) = pair{2};
           else
              error('iosr:resynthesise:unknownParameter','%s is not a recognized parameter name',pair{1})
           end
        end
    end

    delay = options.delay;
    kernel = options.kernel;
    outfilter = options.filter;
    frame_d = options.frame_length;
    order = options.order;

    % validate
    assert(isvector(x) && isnumeric(x), 'iosr:resynthesise:invalidX', 'x must be a vector.');
    assert(all(round(delay)==delay), 'iosr:resynthesise:invalidDelay', 'Values in ''delay'' must be integers.');
    assert(length(delay)==numchans && isnumeric(delay), 'iosr:resynthesise:invalidDelay', 'The ''delay'' parameter must be a numeric vector the same length as the number of columns in m.');
    assert(ischar(outfilter), 'iosr:resynthesise:invalidFilter', 'The ''filter'' option must be a char array (string).');
    assert(round(frame_d)==frame_d && isscalar(frame_d), 'iosr:resynthesise:invalidFrame', '''frame_length'' must be a integer scalar.');
    assert(isnumeric(kernel), 'iosr:resynthesise:invalidKernel', '''kernel'' must be numeric.');
    assert(isvector(cfs) && length(cfs)==numchans && isnumeric(cfs), 'iosr:resynthesise:invalidCfs', 'cfs must be a vector with the same number of elements as there are columns of m.')

    %% get sample-by-sample masks

    if frame_d>1
        m_full = iosr.bss.getFullMask(m,frame_d,delay,kernel);
    else
        m_full = m;
    end

    % crop or zero-pad x to match mask
    if size(m_full,1)<length(x)
        xpad = x(1:length(m_full));
    elseif size(m_full,1)>length(x)
        xpad = [x; zeros(size(m_full,1)-length(x),1)];
    else
        xpad = x;
    end

    %% Pass x through reconstruction filterbank

    switch lower(outfilter)
        case 'gammatone'
            z = iosr.auditory.gammatoneFast(xpad,cfs,fs,true); % phase aligned GTFB
        case 'sinc'
            fcs = [0 cfs];
            z = zeros(length(xpad),numchans);
            for i = 1:numchans
                z(:,i) = iosr.dsp.sincFilter(xpad,fcs(i:i+1)./(fs/2),order)';
            end
        otherwise
            error('Unknown ''filter'' option specified')
    end

    %% create output

    % apply mask
    z = z.*m_full;

    % Sum to create waveforms
    output = sum(z,2);

    % crop or zero-pad output to match input
    if length(output)>length(x)
        output = output(1:length(x));
    elseif length(output)<length(x)
        output = [output; zeros(length(x)-length(output),1)];
    end

end
