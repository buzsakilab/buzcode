function [wavelet, t] = MorletWavelet( f, numcyc, si )
%MorletWavelet(f,numcyc,si) creates a Mortlet Wavelet of numcyc cycles with
%frequency f (Hz).  si: sampling interval
%   Adapted from Cohen Ch13
%
%TO DO:
%   -Add ability to define width in bandwidth instead of cycles (varargin)
%
%Last Updated: 3/15/15
%DLevenstein

% parameters...
s = numcyc/(2*pi*f);    %SD of the gaussian
tbound = (4*s);   %time bounds - at least 4SD on each side, 0 in center
tbound = si*ceil(tbound/si);
t = -tbound:si:tbound;    % time

% and together they make a wavelet
sinusoid = exp(2*pi*1i*f.*t);
gauss = exp(-(t.^2)./(2*s^2));
%A = 1/sqrt(s*sqrt(pi));        %A: normalization... must fix this ASAP
A = 1/(s*sqrt(2*pi));
%see https://groups.google.com/forum/#!topic/analyzingneuraltimeseriesdata/KDJm-Qllq78
%A = ((pi*(2*s^2))^(-0.5));
%A=s^-1 * (2/pi)^(1/2)
%A = sqrt(2*numcyc)*f/(pi*1/si);
A = 1;
%gauss = gauss/(sum(gauss)*s); %normalize to gausswidth.... check this!
wavelet = A * sinusoid .* gauss; 


%normalize to unit energy
wavelet = wavelet./norm(wavelet);

%from matlab
%((pi*FB)^(-0.5))*exp(2*i*pi*FC*X)*exp(-X^2/FB) 

end

