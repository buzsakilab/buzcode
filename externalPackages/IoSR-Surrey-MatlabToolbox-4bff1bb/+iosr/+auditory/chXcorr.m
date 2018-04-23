function [ccg,ic] = chXcorr(hc_L,hc_R,fs,varargin)
%CHXCORR Calculate cross-correlograms with a wide range of options.
% 
%   CCG = IOSR.AUDITORY.CHXCORR(HC_L,HC_R,FS) cross-correlates the input
%   2-D matrices HC_L and HC_R over 10ms frame with a maximum lag of 1ms.
%   It is assumed that the number of frequency channels is min(size(HC_L))
%   and hence HC_L and HC_R can be in either orientation. The
%   cross-correlograms consist of cross-correlations for every frame and
%   frequency channel. CCG has dimensions [lag,frequency,frame]. The
%   function calculates running cross-correlations for every sample and
%   integrates these cross-correlations over each frame. The number of
%   frames frame_count is calculated thus:
% 
%       frame_count = ...
%           floor((max(size(hc_L))-maxlag-1)/frame_length);
% 
%   The underlying cross-correlation algorithm is based on that proposed by
%   Faller & Merimaa [1]. In this implmentation, the time constant of the
%   backward infinite exponential window is given by tau (in samples).
%   
%   CCG = IOSR.AUDITORY.CHXCORR(HC_L,HC_R,FS,'PARAMETER',VALUE) allows a
%   number of options to be specified. The options are:
% 
%   ({} indicates the default value)
% 
%   'frame_length'   : {round(0.01*fs)} | scalar
%       The length of frames used to calculate for integrating
%       cross-correlations.
%   'noverlap'       : {1} | scalar
%       The number of frames over which to integrate the
%       cross-correlations. Note that the frame count is reduced
%       accordingly.
%   'maxlag'         : {round(0.001*fs)} | scalar
%       The maximum lag of the cross-correlation.
%   'tau'            : {round(0.01*fs)} | scalar
%       The time constant of the exponential window used to calculate
%       running cross-correlations. It is recommended that if norm_flag = 1
%       then tau >> 1. As can be seen below, as tau -> 1 then 
%       [aL(m,i,n?1), aR(m,i,n?1)] -> 0, and hence c(m,i,n) -> 1.
%   'inhib'          : {[]} | array
%       Specificies an array with which to multiply the cross-correlations
%       before they are integrated. The value defaults to an empty array,
%       meaning that no inhibition will be applied.
%   'ic_t'           : {0} | scalar
%       Specifies the interaural coherence (IC) threshold. Only samples for
%       which the IC exceeds this threshold will be used to integrate
%       cross-correlations. The algorithm calculates Interaural Coherence
%       (IC) according to [1]. The value should be in the range [0,1].
%   'norm_flag'      : {0} | scalar
%       Specifies whether the cross-correlograms are calculated using
%       normalised cross-correlations. A non-zero value indicates that
%       normalised cross-correlations are used.
%   'inhib_mode'     : {'subtract'} | 'multiply'
%       Specify how the inhibition is applied. The default 'subtract' will
%       subtract inhib from the running cross-correlations (negative values
%       are set to zero); 'multiply' will multiply inhib with the running
%       cross-correlations.
% 
%   [CCG,IC] = IOSR.AUDITORY.CHXCORR(...) returns the calculated IC to the
%   matrix IC. Although the matrix returned is the same size as hc_L, IC is
%   only calculated for samples 1:frame_count*frame_length, other values
%   will be set to 0.
% 
%   Algorithm
% 
%   The running normalised cross-correlation is calculated as [1]:
% 
%   C(m,i,n) = c(m,i,n) / sqrt( aL(m,i,n) * aR(m,i,n)
% 
%   where
% 
%   c(m,i,n) = (1/tau) * HC_L(i, max(n+m,n)) * HC_R(i, max(n-m,n)) + ...
%       (1-1/tau) * c(m,i,n-1),
% 
%   aL(m,i,n) = (1/tau) * (HC_L(i, max(n+m,n)))^2 + ...
%       (1-1/tau) * aL(m,i,n-1),
% 
%   aR(m,i,n) = (1/tau) * (HC_R(i, max(n+m,n)))^2 + ...
%       (1-1/tau) * aR(m,i,n-1),
% 
%   i is the frequency index, m is the lag index, and n is the sample
%   index. The running (non-normalised) cross-correlation is calculated is
%   simply c(m,i,n).
% 
%   The interaural coherence is
% 
%   IC(i,n) = max(C(m,i,n)) % (i.e. over m)
% 
%   The cross-correlogram is calculated as the sum of cross- correlations
%   in a given frame:
% 
%   CCG(m,i,j) = sum(C(m,i,J),3)
% 
%   where
% 
%   J = (j-1) * frame_length + 1 : j * frame_length
% 
%   References
% 
%   [1] C. Faller and J. Merimaa, "Source localization in complex listening
%       situations: Selection of binaural cues based on interaural
%       coherence", The Journal of the Acoustical Society of America, vol.
%       116, pp.3075-3089, Nov. 2004.
% 
%   Further Reading
%   
%   C. Hummersone, R. Mason, and T. Brookes, "A comparison of computational
%       precedence models for source separation in reverberant
%       environments", The Journal of the Audio Engineering Society, vol.
%       61(7/8), pp.508-520, July 2013.

%   Copyright 2016 University of Surrey.

    assert(nargin>=3, 'iosr:chXcorr:nargin', 'Number of input arguments must be greater than or equal to three.')

    if isparameter(varargin,'inhib_mode') && ~isparameter(varargin,'inhib')
        warning('iosr:instIld:inhibMode','''inhib_mode'' specified, but no inhibition array ''inhib''.');
    end

    % Check source file is compiled
    iosr.general.checkMexCompiled('-largeArrayDims',fullfile(fileparts(mfilename('fullpath')),'chXcorr_c.c'))

    options = struct(...
        'frame_length',round(0.01*fs),...
        'noverlap',1,...
        'maxlag',round(0.001*fs),...
        'tau',round(0.01*fs),...
        'inhib',[],...
        'ic_t',0,...
        'norm_flag',0,...
        'inhib_mode','subtract');

    % read parameter/value inputs
    if nargin > 3 % if parameters are specified
        % read the acceptable names
        optionNames = fieldnames(options);
        % count arguments
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('iosr:chXcorr:nameValuePair','CHXCORR needs propertyName/propertyValue pairs')
        end
        % overwrite defults
        for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
           IX = strcmpi(pair{1},optionNames); % find match parameter names
           if any(IX)
              % do the overwrite
              options.(optionNames{IX}) = pair{2};
           else
              error('iosr:chXcorr:unknownOption','%s is not a recognized parameter name',pair{1})
           end
        end
    end

    % assign options to variables
    frame_length = options.frame_length;
    noverlap = options.noverlap;
    maxlag = options.maxlag;
    tau = options.tau;
    inhib_mode = options.inhib_mode;
    norm_flag = options.norm_flag;
    ic_t = options.ic_t;
    inhib = options.inhib;

    % check inputs
    assert(all(size(hc_L)==size(hc_R)), 'iosr:chXcorr:invalidInput', '''hc_L'' and ''hc_R'' must be the same size')
    assert(round(frame_length)==frame_length && isscalar(frame_length) && frame_length>0, 'iosr:chXcorr:invalidFrame', ...
        '''frame_length'' must be an integer greater than zero')
    assert(round(noverlap)==noverlap && isscalar(noverlap) && noverlap>0, 'iosr:chXcorr:invalidNoverlap', ...
        '''noverlap'' must be an integer greater than zero')
    assert(round(maxlag)==maxlag && isscalar(maxlag) && maxlag>0, 'iosr:chXcorr:invalidMaxlag', ...
        '''maxlag'' must be an integer greater than zero')
    assert(isscalar(tau) && tau>=1, 'iosr:chXcorr:invalidTau', '''tau'' must be a scalar greater than or equal to one')
    assert(isscalar(norm_flag), 'iosr:chXcorr:invalidNorm', '''norm_flag'' must be a scalar')
    assert(isscalar(ic_t) && ic_t>=0 && ic_t<=1, 'iosr:chXcorr:invalidIct', '''ic_t'' must be a scalar in the range [0,1]')
    assert(ischar(inhib_mode), 'iosr:chXcorr:invalidInhibMode', '''inhib_mode'' must be a char array (string)')

    % Calculate frame count
    frame_count = floor(max(size(hc_L))/(frame_length));
    frame_count = frame_count-noverlap+1;

    % Calculate number of frequency channels
    numchans = min(size(hc_L));
    numsamples = max(size(hc_L));

    % Check orientation of HC and inhib data (i.e. that frequency runs across the rows)
    dims = size(hc_L);
    hc_L = check_input(hc_L,2,numchans);
    hc_R = check_input(hc_R,2,numchans);

    % set a flag if data has been transposed in this way
    if dims(1)~=size(hc_L,1)
        rot = true;
    else
        rot = false;
    end

    % set inhibition mode ID
    switch inhib_mode
        case 'multiply'
            inhib_mode_ID = 1;
            if isempty(inhib)
                inhib = ones(size(hc_L));
            end
        case 'subtract'
            inhib_mode_ID = 2;
            if isempty(inhib)
                inhib = zeros(size(hc_L));
            end
        otherwise
            error('iosr:chXcorr:unknownInhibMode','''inhib_mode'' must be set to ''multiply'' or ''subtract''')
    end

    inhib = check_input(inhib,2,numchans);

    % Append HC and inhibition data with zeros for cross-correlation
    hc_L = [hc_L; zeros(maxlag+1,numchans)];
    hc_R = [hc_R; zeros(maxlag+1,numchans)];
    inhib = [inhib; zeros(maxlag+1,numchans)];

    assert(all(size(inhib)==size(hc_L)), 'iosr:chXcorr:invalidInhib', '''inhib'' must be a matrix the same size as ''hc_L'' or ''hc_R''')

    % Calculate cross-correlograms
    [ccg,ic] = iosr.auditory.chXcorr_c(hc_L,hc_R,frame_count,frame_length,noverlap,maxlag,tau,inhib,ic_t,norm_flag,inhib_mode_ID);

    % Correct orientation of IC data, if data was transposed, and crop to remove appended zeros
    ic = ic(1:numsamples,:);
    if rot
        ic = ic';
    end

end

function output = check_input(input,dim,target)
%CHECK_INPUT check input is correct orientation

    if size(input,dim)~=target
        output = input';
        assert(size(output,dim)==target, 'iosr:chXcorr:invalidInputs', 'Input invalid')
    else
        output = input;
    end

end

function set = isparameter(input,parameter)
%ISPARAMETER check for input parameter

    set = any(strcmpi(input(cellfun(@ischar,input)),parameter));

end
