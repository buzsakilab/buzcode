function [ convolution_result_fft ] = FConv( kernel, signal, varargin )
%FConv(kernel,signal) convolves a signal with a kernal via the fourier
%transform for speedy delivery.
%Note: as is, kernel is centered. 
%   Adapted from Cohen Chapter 13
%
%TO DO:
%   -add varargin in for centered or left-aligned kernel
%
%Last Updated: 7/30/15
%DLevenstein

% Kernel and signal must be same orientation


% FFT parameters
n_kernel            = length(kernel);
n_signal               = length(signal);
n_convolution        = n_kernel + n_signal-1;
half_of_kernel_size = (n_kernel-1)/2;


% FFT of wavelet and EEG data
fft_kernel = fft(kernel,n_convolution);
fft_signal = fft(signal,n_convolution);

convolution_result_fft = ifft((fft_kernel.*fft_signal),n_convolution);

% cut off edges
convolution_result_fft = convolution_result_fft(half_of_kernel_size+1:end-half_of_kernel_size);


end