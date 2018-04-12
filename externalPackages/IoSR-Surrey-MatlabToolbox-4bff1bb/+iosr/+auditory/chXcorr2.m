function [ccg,ic] = chXcorr2(hc_L,hc_R,fs,varargin)
%CHXCORR2 Calculate cross-correlograms with a range of options.
% 
%   CCG = IOSR.AUDITORY.CHXCORR2(HC_L,HC_R,FS) cross-correlates the input
%   2-D matrices HC_L and HC_R over 10ms frame with a maximum lag of 1ms.
%   It is assumed that the number of frequency channels is min(size(HC_L))
%   and hence HC_L and HC_R can be in either orientation. The
%   cross-correlograms consist of cross-correlations for every frame and
%   frequency channel. CCG has dimensions [lag,frequency,frame]. The
%   function calculates the traditional cross-correlation in each frame.
%   The number of frames FRAME_COUNT is calculated thus:
% 
%       FRAME_COUNT = FIX((MAX(SIZE(HC_L)))/FRAME_LENGTH);
%   
%   CCG = IOSR.AUDITORY.CHXCORR2(HC_L,HC_R,FS,'PARAMETER',VALUE) allows a
%   number of options to be specified. The options are:
% 
%   ({} indicates the default value)
% 
%   'frame_length'   : {round(0.01*fs)} | scalar
%       The length of frames (in samples) used for calculating
%       cross-correlations.
%   'hop'            : {[]} | scalar
%       The hop size (in samples). By default the hop size is equal to the
%       frame length (HOP is specified as an empty array). The hop size
%       determines the number of frames as
%           FIX((MAX(SIZE(HC_L))-(FRAME_LENGTH-HOP))/HOP);
%   'maxlag'         : {round(0.001*fs)} | scalar
%       The maximum lag of the cross-correlation (in samples).
%   'norm_flag'      : {0} | scalar
%       Specifies whether the cross-correlograms are calculated using
%       normalised cross-correlations. A non-zero value indicates that
%       normalised cross-correlations are used.
% 
%   [CCG,IC] = IOSR.AUDITORY.CHXCORR2(...) returns the calculated IC for
%   each frame to the matrix IC.

%   Copyright 2016 University of Surrey.

    assert(nargin>=3, 'iosr:chXcorr2:nargin', 'Number of input arguments must be greater than or equal to three.')

    % Check source file is compiled
    iosr.general.checkMexCompiled('-largeArrayDims',fullfile(fileparts(mfilename('fullpath')),'chXcorr2_c.c'))

    options = struct(...
        'frame_length',round(0.01*fs),...
        'maxlag',round(0.001*fs),...
        'norm_flag',0,...
        'hop',[]);

    % read parameter/value inputs
    if nargin > 3 % if parameters are specified
        % read the acceptable names
        optionNames = fieldnames(options);
        % count arguments
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('iosr:chXcorr2:nameValuePair','CHXCORR2 needs propertyName/propertyValue pairs')
        end
        % overwrite defults
        for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
           IX = strcmpi(pair{1},optionNames); % find match parameter names
           if any(IX)
              % do the overwrite
              options.(optionNames{IX}) = pair{2};
           else
              error('iosr:chXcorr2:unknownOption','%s is not a recognized parameter name',pair{1})
           end
        end
    end

    % assign options to variables
    frame_length = options.frame_length;
    maxlag = options.maxlag;
    norm_flag = options.norm_flag;
    if isempty(options.hop)
        hop = frame_length;
    else
        hop = options.hop;
    end

    % check inputs
    assert(all(size(hc_L)==size(hc_R)), 'iosr:chXcorr2:invalidInput', '''hc_L'' and ''hc_R'' must be the same size')
    assert(round(frame_length)==frame_length && isscalar(frame_length) && frame_length>0, ...
        'iosr:chXcorr2:invalidFrame', '''frame_length'' must be an integer greater than zero')
    assert(round(maxlag)==maxlag && isscalar(maxlag) && maxlag>0, 'iosr:chXcorr2:invalidMaxlag', ...
        '''maxlag'' must be an integer greater than zero')
    assert(isscalar(norm_flag), 'iosr:chXcorr2:invalidNorm', '''norm_flag'' must be a scalar')

    % Calculate frame count
    frame_count = fix((max(size(hc_L))-(frame_length-hop))/hop);

    % Calculate number of frequency channels
    numchans = min(size(hc_L));

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

    % Calculate cross-correlograms
    [ccg,ic] = iosr.auditory.chXcorr2_c(hc_L,hc_R,frame_count,frame_length,maxlag,hop,norm_flag);

    % Correct orientation of IC data, if data was transposed, and crop to remove appended zeros
    if rot
        ic = ic';
    end

end

function output = check_input(input,dim,target)
%CHECK_INPUT check input is correct orientation

    if size(input,dim)~=target
        output = input';
        assert(size(output,dim)==target, 'iosr:chXcorr2:invalidInput', 'Input invalid')
    else
        output = input;
    end

end
