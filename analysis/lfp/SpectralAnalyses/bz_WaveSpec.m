function [freqs,t,spec] = bz_WaveSpec(signal,frange,nfreqs,ncyc,si,space)
%[freqs,t,spec] = WaveSpec(signal,frange,nfreqs,ncyc,si,space) calculates the 
%wavelet transform of a signal with nfreqs frequencies in the range frange 
%[fmin fmax]. Spacing between frequencies can be 'lin' or 'log'.
%Time-frequency resolution is defined by ncyc, the number of cycles in each
%wavelet. Uses Morlet (Gabor) wavelet.
%
%
%INPUT
%   signal  [Nt x 1] 
%       (Optional) cell array of [Nt x 1] time series
%   frange  [low frequency, high frequency]
%   nfreqs  number of frequencies to divide your signal into
%   ncyc    number of cycles in your wavelet (recommend: 5)
%   si      sampling interval (1/sampling frequency)
%   space   'log' or 'lin'  spacing of your frequencies
%
%OUTPUT
%   freqs
%   t
%   spec    spectrogram!
%
%
%
%
%TO DO:
%   -Log frequency response wavelets? (new function)
%   -Matricise For loop?
%   -Add LaTeX formula to comments
%   -Don't need to FFT every freq...
%   -update t output for if LFP is a cell array
%
%
%Dependencies
%   WaveFilt
%   MorletWavelet
%   FConv
%
%
%Last Updated: 10/9/15
%DLevenstein
%%
if isempty(signal)
    freqs=[]; t=[]; spec=signal;
    return
end

if iscell(signal)
    celllengths = cellfun(@length,signal);
    signal = vertcat(signal{:});
end

if ~isa(signal,'single') | ~isa(signal,'double')
    signal = single(signal);
end


fmin = frange(1);
fmax = frange(2);
if strcmp(space,'log')
    assert(fmin~=0,'Log-spaced frequencies cannot have min of 0')
    freqs = logspace(log10(fmin),log10(fmax),nfreqs);
elseif strcmp(space,'lin')
    freqs = linspace(fmin,fmax,nfreqs);
else
    display('Frequency spacing must be "lin" or "log".')
end

t = (0:length(signal)-1)*si;

spec = zeros(nfreqs,length(t));
for f_i = 1:nfreqs
    if mod(f_i,10) == 1;
        display(['freq ',num2str(f_i),' of ',num2str(nfreqs)]);
    end
    spec(f_i,:) = WaveFilt(signal,freqs(f_i),ncyc,si);
end

if exist('celllengths','var')
    spec = mat2cell(spec,nfreqs,celllengths);
end
    

end

