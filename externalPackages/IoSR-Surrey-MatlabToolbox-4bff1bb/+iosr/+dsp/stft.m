function [s,f,t] = stft(x,nfft,hop,fs)
%STFT Calculate the short-time Fourier transform of a signal
% 
%   S = IOSR.DSP.STFT(X) calculates the short-time Fourier transform (STFT)
%   of signal X using 1024-point segments in steps (hops) of 512 points. S
%   is a series of one-sided FFTs of size [N,K] where N is the number of
%   frequency bins and K is the number of time frames, such that N=513 and
%   K=FIX((LENGTH(X)-512)/512). A hamming window is applied to each time
%   segment. The window is normalised such that it fulfills the constant
%   overlap-add (COLA) criterion. Fulfilling the criterion is useful when
%   resynthesising the signal using the ISTFT, as the magnitude of the
%   output of the ISTFT will be retained.
% 
%   S = IOSR.DSP.STFT(X,NFFT) calculates the STFT using NFFT-point
%   segments. The hop size is HOP=FIX(NFFT/2);
%   K=FIX((LENGTH(X)-(NFFT-HOP))/HOP), and N depends on whether NFFT is odd
%   or even. If NFFT is even, N=(NFFT/2)+1; if NFFT is odd, N=(NFFT+1)/2.
% 
%   S = IOSR.DSP.STFT(X,WINDOW) calculates the STFT with
%   NFFT=LENGTH(WINDOW) and applies the vector WINDOW to each segment. The
%   window is normalised such that it fulfills the constant overlap-add
%   (COLA) criterion. The normalised window is calculated as
%   HOP.*WINDOW./SUM(WINDOW).
% 
%   S = IOSR.DSP.STFT(X,WINDOW,HOP) steps through X by HOP samples.
% 
%   [S,W] = IOSR.DSP.STFT(...) returns the normalised frequencies for each
%   bin to W.
% 
%   [S,F] = IOSR.DSP.STFT(X,WINDOW,HOP,FS) returns the frequencies for each
%   bin to F according to sampling frequency FS.
% 
%   [S,F,T] = IOSR.DSP.STFT(...) returns the corresponding time T (seconds)
%   for each column of S.
% 
%   Example
% 
%       % plot a spectrogram
%       load handel.mat
%       [Y,f,t] = iosr.dsp.stft(x,1024,128,Fs);
%       figure
%       imagesc(t,f,20.*(log10(abs(Y))));
%       set(gca,'ydir','normal')
%       ylabel('Frequency [Hz]')
%       xlabel('Time [s]')
% 
%   See also IOSR.DSP.ISTFT.

%   Copyright 2016 University of Surrey.

    %% check input
    
    assert(isvector(x) && numel(x)>1, 'iosr:stft:invalidX', 'X must be a vector')
    
    % check nfft
    if nargin<2
        nfft = 1024;
    end
    
    % determine window
    if numel(nfft)>1
        win = nfft;
        assert(isvector(win), 'iosr:stft:invalidNfft', 'WINDOW must be a vector')
        nfft = length(win);
    else
        assert(round(nfft)==nfft && nfft>0, 'iosr:stft:invalidNfft', 'NFFT must be a positive integer')
        win = hamming(nfft);
    end
    
    % check x length
    assert(length(x)>=nfft, 'iosr:stft:invalidInput', 'X must have at least NFFT samples')
    
    % determine hop
    if nargin<3
        hop = fix(nfft/2);
    else
        assert(isscalar(hop) & round(hop)==hop, 'iosr:stft:invalidHop', 'HOP must be an integer')
        assert(hop<=nfft && hop>0, 'iosr:stft:invalidHop', 'HOP must be less than or equal to NFFT, and greater than 0')
    end
    
    % normalise window
    win = hop.*win./sum(win);
    
    % determine fs
    if nargin<4
        fs = 1;
    else
        assert(isscalar(fs), 'iosr:stft:invaldFs', 'FS must be an scalar')
    end
    
    %% calculate outputs
    
    % calculate highest bin for one-sided FFT
    if mod(nfft,2)==0
        Nout = (nfft/2)+1;
    else
        Nout = (nfft+1)/2;
    end
    
    % calculate number of frames
    K = fix((length(x)-(nfft-hop))/hop);
    
    % calculate STFTs
    s = zeros(Nout,K);
    for k = 1:K
        samples = ((k-1)*hop)+1:((k-1)*hop)+nfft;
        temp = fft(win.*x(samples));
        s(:,k) = temp(1:Nout);
    end
    
    % calculate frequency and time
    if nargout>1
        f = fs.*((0:Nout-1)./nfft)';
    end
    if nargout>2
        t = (0:K-1).*(hop/fs);
    end

end
