function [ filt ] = WaveFilt( signal, f, numcyc, si )
%WaveFilt(signal,f,numcyc,si) filters signal with a mortlet wavelet of
%numcyc cycles of frequency f.
%
%TO DO:
%   -Add option to define wavelet width in bandwidth instead of cycles
%   -Optional Arguement: 'fourier' - signal is already fourier transformed
%       for WaveSpec
%
%Last Updated: 4/5/15
%DLevenstein

wavelet = MorletWavelet(f,numcyc,si);
filt = FConv(wavelet',signal);


end

