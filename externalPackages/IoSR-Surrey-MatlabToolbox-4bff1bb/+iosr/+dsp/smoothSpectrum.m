function x_oct = smoothSpectrum(X,f,Noct)
%SMOOTHSPECTRUM Apply 1/N-octave smoothing to a frequency spectrum
% 
%   X_OCT = IOSR.DSP.SMOOTHSPECTRUM(X,F,NOCT) applies 1/NOCT-octave
%   smoothing to the frequency spectrum contained in vector X sampled at
%   frequencies in vector F. X can be a log-, magnitude-, or
%   power-spectrum. Setting Noct to 0 results in no smoothing.
%   
%   Algorithm
%   
%   The function calculates the i-th smoothed spectral coefficient X_OCT(i)
%   as the sum of the windowed spectrum. The window is a Gaussian whose
%   centre frequency is F(i), and whose standard deviation is proportional
%   to F(i)/NOCT.
% 
%   Example
% 
%       % Calculate the 1/3-octave-smoothed power spectral density of the
%       % Handel example.
% 
%       % load signal
%       load handel.mat
%       
%       % take fft
%       Y = fft(y);
%       
%       % keep only meaningful frequencies
%       NFFT = length(y);
%       if mod(NFFT,2)==0
%           Nout = (NFFT/2)+1;
%       else
%           Nout = (NFFT+1)/2;
%       end
%       Y = Y(1:Nout);
%       f = ((0:Nout-1)'./NFFT).*Fs;
%       
%       % put into dB
%       Y = 20*log10(abs(Y)./NFFT);
%       
%       % smooth
%       Noct = 3;
%       Z = iosr.dsp.smoothSpectrum(Y,f,Noct);
%       
%       % plot
%       figure
%       semilogx(f,Y,f,Z)
%       grid on
% 
%   See also IOSR.DSP.LTAS, FFT.

%   Copyright 2016 University of Surrey.

    %% Input checking

    assert(isvector(X), 'iosr:smoothSpectrum:invalidX', 'X must be a vector.');
    assert(isvector(f), 'iosr:smoothSpectrum:invalidF', 'F must be a vector.');
    assert(isscalar(Noct), 'iosr:smoothSpectrum:invalidNoct', 'NOCT must be a scalar.');
    assert(isreal(X), 'iosr:smoothSpectrum:invalidX', 'X must be real.');
    assert(all(f>=0), 'iosr:smoothSpectrum:invalidF', 'F must contain positive values.');
    assert(Noct>=0, 'iosr:smoothSpectrum:invalidNoct', 'NOCT must be greater than or equal to 0.');
    assert(isequal(size(X),size(f)), 'iosr:smoothSpectrum:invalidInput', 'X and F must be the same size.');

    %% Smoothing
    
    % calculates a Gaussian function for each frequency, deriving a
    % bandwidth for that frequency

    x_oct = X; % initial spectrum
    if Noct > 0 % don't bother if no smoothing
        for i = find(f>0,1,'first'):length(f)
            g = gauss_f(f,f(i),Noct);
            x_oct(i) = sum(g.*X); % calculate smoothed spectral coefficient
        end
        % remove undershoot when X is positive
        if all(X>=0)
            x_oct(x_oct<0) = 0;
        end
    end

end

function g = gauss_f(f_x,F,Noct)
% GAUSS_F calculate frequency-domain Gaussian with unity gain
% 
%   G = GAUSS_F(F_X,F,NOCT) calculates a frequency-domain Gaussian function
%   for frequencies F_X, with centre frequency F and bandwidth F/NOCT.

    sigma = (F/Noct)/pi; % standard deviation
    g = exp(-(((f_x-F).^2)./(2.*(sigma^2)))); % Gaussian
    g = g./sum(g); % normalise magnitude

end
