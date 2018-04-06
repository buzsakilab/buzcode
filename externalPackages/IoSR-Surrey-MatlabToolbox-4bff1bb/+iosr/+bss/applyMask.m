function [z,t] = applyMask(s,m,nfft,hop,fs)
%APPLYMASK Apply a time-frequency mask to an STFT
% 
%   Z = IOSR.BSS.APPLYMASK(S,M) applies the time-frequency mask M to STFT
%   S. S has dimensions [N,K] where N is the number of frequency bins and K
%   is the number of time frames. The ISTFT is calculated using 1024-point
%   FFTs and hop sizes of 512 points.
% 
%   Z = IOSR.BSS.APPLYMASK(S,M,NFFT) uses NFFT-length segments in the
%   ISTFT.
% 
%   Z = IOSR.BSS.APPLYMASK(S,M,NFFT,HOP) uses hop size HOP for the ISTFT.
% 
%   [Z,T] = IOSR.BSS.APPLYMASK(S,M,NFFT,HOP,FS) uses sampling frequency FS
%   to return the corresponding time T of each element in Z.
% 
%   See also IOSR.DSP.STFT, IOSR.DSP.ISTFT, IOSR.BSS.IDEALMASKS, 
%            IOSR.BSS.APPLYIDEALMASKS.

%   Copyright 2016 University of Surrey.

    %% check input
    
    % check input
    if nargin<2
        error('iosr:applyMask:nargin','Not enough input arguments')
    end
    
    % check compulsory inputs
    assert(isequal(size(s),size(m)), 'iosr:applyMask:invalidInput', 'S and M must be the same size')
    
    % check nfft
    if nargin<3
        nfft = 1024;
    else
        assert(isscalar(nfft) && round(nfft)==nfft && nfft>0, 'iosr:applyMask:invalidNfft', 'NFFT must be a positive scalar integer')
    end
    
    % determine hop
    if nargin<4
        hop = fix(nfft/2);
    else
        assert(isscalar(hop) & round(hop)==hop, 'iosr:applyMask:invalidHop', 'HOP must be an integer')
        assert(hop<=nfft && hop>0, 'iosr:applyMask:invalidHop', 'HOP must be less than or equal to NFFT, and greater than 0')
    end
    
    % determine fs
    if nargin<5
        fs = 1;
    else
        assert(isscalar(fs), 'iosr:applyMask:invalidFs', 'FS must be an scalar')
    end
    
    %% calculate outputs

    % apply mask and return signal
    z = iosr.dsp.istft(s.*m,nfft,hop);
    
    % return time
    if nargout>1
        t = (0:length(z)-1)./fs;
    end
    
end
