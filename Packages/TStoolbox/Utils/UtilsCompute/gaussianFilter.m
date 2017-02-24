function y = gaussianFilter(x,fs,fc)

% y = gaussianFilter(x,fs,fc)
% INPUTS:
%   x: vector (or column matrix) of data
%   fs: sampling fq
%   fc: high pass cut off frequency or a 2 value vector of [lowpass
%   highpass] frequencies of the passband filter.

% drein Peyrache

sigma = fs./(2*pi*fc);
y = gaussFilt(x,sigma(1),0);

if length(fc)==2
    y = gaussFilt(x,sigma(2),0)-y;
elseif length(fc)>2
    error('Cut-off frequency has either one or two values')
end
