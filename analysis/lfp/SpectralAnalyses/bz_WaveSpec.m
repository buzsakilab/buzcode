function [wavespec] = bz_WaveSpec(lfp,varargin)
%[wavespec] = bz_WaveSpec(lfp) calculates the 
%wavelet transform of a signal with nfreqs frequencies in the range frange 
%[fmin fmax]. Spacing between frequencies can be 'lin' or 'log'.
%Time-frequency resolution is defined by ncyc, the number of cycles in each
%wavelet. Uses Morlet (Gabor) wavelet.
%
%
%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       'frange'	[low frequency, high frequency]     (default: [1 128])
%       'nfreqs'  	number of frequencies               (default: 100
%       'ncyc'      number of cycles in your wavelet    (default: 5)
%       'space'     'log' or 'lin'  spacing of f's      (default: 'log')
%       'samplingRate' (only if input is not a buzcode structure)
%       'intervals'  ADD THIS - ability to spec intervals
%       'showprogress' true/false (default:false)
%       'saveMat'   put the basePath to save an LFP file
%    =========================================================================
%
%OUTPUT
%   wavespec            buzcode-style structure
%       .data           [t x nfreqs] your spectrogram
%       .timestamps     [t x 1] timestamps
%       .freqs          frequencies of each column
%       .samplingRate   (Hz)
%       .filterparms    a structure that holds the parameters used for
%                       filtering, for future reference
%

%TO DO:
%   -Matricise For loop?
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
%% Parse the inputs

%Parameters
parms = inputParser;
addParameter(parms,'frange',[1 128],@isnumeric);
addParameter(parms,'nfreqs',100,@isnumeric);
addParameter(parms,'ncyc',5,@isnumeric);
addParameter(parms,'space','log');
addParameter(parms,'samplingRate',[]);
addParameter(parms,'showprogress',false,@islogical);
addParameter(parms,'saveMat',false);

parse(parms,varargin{:})
frange = parms.Results.frange;
nfreqs = parms.Results.nfreqs;
ncyc = parms.Results.ncyc;
space = parms.Results.space;
samplingRate = parms.Results.samplingRate;
showprogress = parms.Results.showprogress;
saveMat = parms.Results.saveMat;


%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif isempty(lfp)
    wavespec = lfp;
    return
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

si = 1./samplingRate;


%%




if ~isa(data,'single') || ~isa(data,'double')
    data = single(data);
end

%Frequencies
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


%Filter with wavelets
spec = zeros(length(timestamps),nfreqs);
for f_i = 1:nfreqs
    if showprogress
        if mod(f_i,10) == 1
            display(['freq ',num2str(f_i),' of ',num2str(nfreqs)]);
        end  
    end
    wavelet = MorletWavelet(freqs(f_i),ncyc,si);
    spec(:,f_i) = FConv(wavelet',data);
end

if exist('celllengths','var')
    spec = mat2cell(spec,nfreqs,celllengths);
end

%% Output in buzcode format

wavespec.data = spec;
wavespec.timestamps = timestamps;
wavespec.freqs = freqs;
wavespec.samplingRate = samplingRate;
wavespec.filterparms.ncyc = ncyc;
wavespec.filterparms.nfreqs = nfreqs;
wavespec.filterparms.frange = frange;
wavespec.filterparms.space = space;

if saveMat
    baseName = bz_BasenameFromBasepath(saveMat);
    lfpfilename = fullfile(basePath,[baseName,'.wavespec.lfp.mat']);
    save(lfpfilename,wavespec)
end

end

