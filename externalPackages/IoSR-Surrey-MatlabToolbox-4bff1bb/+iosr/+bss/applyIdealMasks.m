function [z_irm,z_ibm,t] = applyIdealMasks(xt,xi,nfft,hop,fs)
%APPLYIDEALMASKS Calculate and apply ideal masks via STFT
% 
%   Z_IRM = IOSR.BSS.APPLYIDEALMASKS(XT,XI) calculates the ideal ratio mask
%   (IRM) and applies it to the mixture XT+XI, where XT is the target
%   signal and XI is the interference signal. The IRM is calculated via the
%   STFT using 1024-point windows with 512-point overlap. Z_IRM, XT, and XI
%   are vectors. If XT and XI are of different lengths then the shorter
%   signal is zero-padded in order to make them the same length; Z_IRM is
%   the same length as XT and XI.
% 
%   Z_IRM = IOSR.BSS.APPLYIDEALMASKS(XT,XI,NFFT) uses NFFT-length segments
%   in the STFT.
% 
%   Z_IRM = IOSR.BSS.APPLYIDEALMASKS(XT,XI,WINDOW) uses
%   LENGTH(WINDOW)-length segments in the STFT and applies WINDOW to each
%   segment.
% 
%   Z_IRM = IOSR.BSS.APPLYIDEALMASKS(XT,XI,WINDOW,HOP) uses hop size HOP
%   for the STFT.
% 
%   [Z_IRM,Z_IBM] = IOSR.BSS.APPLYIDEALMASKS(...) calculates the ideal
%   binary mask (IBM) and applies it to the mixture, returning the result
%   to Z_IBM.
% 
%   [Z_IRM,Z_IBM,T] = IOSR.BSS.APPLYIDEALMASKS(XT,XI,WINDOW,HOP,FS) uses
%   sampling frequency FS to return the corresponding time T of each
%   element in Z_IRM and Z_IBM.
% 
%   See also IOSR.DSP.STFT, IOSR.DSP.ISTFT, IOSR.BSS.IDEALMASKS,
%            IOSR.BSS.APPLYMASKS.

%   Copyright 2016 University of Surrey.
    
    %% check input
    
    % check signals
    assert(isvector(xt) && numel(xt)>1, 'iosr:applyIdealMasks:invalidXt', 'XT must be a vector')
    assert(isvector(xi) && numel(xi)>1, 'iosr:applyIdealMasks:invalidXi', 'XI must be a vector')
    
    % make equal length
    maxlength = max([length(xi) length(xt)]);
    xt = pad(xt,maxlength);
    xi = pad(xi,maxlength);
    
    % check nfft
    if nargin<3
        nfft = 1024;
    end
    
    % determine window
    if numel(nfft)>1
        win = nfft;
        assert(isvector(win), 'iosr:applyIdealMasks:invalidWin', 'WINDOW must be a vector')
        nfft = length(win);
    else
        assert(round(nfft)==nfft && nfft>0, 'iosr:applyIdealMasks:invalidNfft', 'NFFT must be a positive integer')
        win = hamming(nfft);
    end
    
    % check x length
    assert(length(xt)>=nfft, 'iosr:applyIdealMasks:invalidXt', 'XT must have at least NFFT samples')
    assert(length(xi)>=nfft, 'iosr:applyIdealMasks:invalidXi', 'XI must have at least NFFT samples')
    
    % determine hop
    if nargin<4
        hop = fix(nfft/2);
    else
        assert(isscalar(hop) & round(hop)==hop, 'iosr:applyIdealMasks:invalidHop', 'HOP must be an integer')
        assert(hop<=nfft && hop>0, 'iosr:applyIdealMasks:invalidHop', 'HOP must be less than or equal to NFFT, and greater than 0')
    end
    
    % determine fs
    if nargin<5
        fs = 1;
    else
        assert(isscalar(fs), 'iosr:applyIdealMasks:invalidFs', 'FS must be an scalar')
    end
    
    %% calculate outputs
    
    % STFTs of signals and mixture
    st = iosr.dsp.stft(xt,win,hop);
    si = iosr.dsp.stft(xi,win,hop);
    mix = iosr.dsp.stft(xt+xi,win,hop);
    
    % return ideal masks
    [irm,ibm] = iosr.bss.idealMasks(st,si);
    
    % apply IRM
    z_irm = iosr.bss.applyMask(mix,irm,nfft,hop,fs);
    z_irm = pad(z_irm,maxlength);
    
    % apply IBM
    if nargout>1
        z_ibm = iosr.bss.applyMask(mix,ibm,nfft,hop,fs);
        z_ibm = pad(z_ibm,maxlength);
    end
    
    % calculate t
    if nargout>2
        t = (0:length(z_irm)-1)./fs;
    end
    
end

function y = pad(x,dur)
%PAD Zero-pad a vector

    if length(x)<dur
        y = [x(:); zeros(dur-length(x),1)];
    else
        y = x(:);
    end

end
