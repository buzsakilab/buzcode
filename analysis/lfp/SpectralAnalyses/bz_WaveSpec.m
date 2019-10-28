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
%       'frange'    [low frequency, high frequency]     (default: [1 128])
%       'nfreqs'    number of frequencies               (default: 100
%       'roundfreqs' round freqs to unique integer vals (default: false
%                       *Note this may decrease number of freqs
%       'nfreqs'    number of frequencies               (default: 100
%       'ncyc'      number of cycles in your wavelet    (default: 5)
%       'fvector'   predefined vector of frequencies 
%       'space'     'log' or 'lin'  spacing of f's      (default: 'log')
%       'samplingRate' (only if input is not a buzcode structure)
%       'intervals'  restrict your spectrogram to timestamps in specific
%                   intervals
%       'chanID'    if lfp structure has multiple channels, which one would
%                   you like to calcaulte the wavelet transform of?
%                   (note: requires field lfp.channels)
%       'showprogress' true/false (default:false)
%       'saveMat '   put the basePath to save an LFP file
%       'MatNameExtraText'   text(X) to add to name as in: 'basename.wavespec(text).lfp.mat'
%       'downsampleout' factor by which to downsample output (default: 1)
%    =========================================================================
%
%OUTPUT
%   wavespec            buzcode-style structure
%       .data           [t x nfreqs] your spectrogram
%       .timestamps     [t x 1] timestamps
%       .freqs          frequencies of each column
%       .samplingRate   (Hz)
%       .channels       Channel indices of channels filtered... taken from lfp input
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
%Modified by Antonio FR, 7/18/18

%% Parse the inputs

%Parameters
parms = inputParser;
addParameter(parms,'frange',[1 128],@isnumeric);
addParameter(parms,'nfreqs',100,@isnumeric);
addParameter(parms,'ncyc',5,@isnumeric);
addParameter(parms,'space','log');
addParameter(parms,'samplingRate',[]);
addParameter(parms,'showprogress',false,@islogical);
addParameter(parms,'roundfreqs',false,@islogical);
addParameter(parms,'saveMatPath',[]);
addParameter(parms,'MatNameExtraText',[]);
addParameter(parms,'fvector',[]);
addParameter(parms,'intervals',[-Inf Inf])
addParameter(parms,'chanID',[])
addParameter(parms,'downsampleout',1)

parse(parms,varargin{:})
frange = parms.Results.frange;
nfreqs = parms.Results.nfreqs;
ncyc = parms.Results.ncyc;
space = parms.Results.space;
samplingRate = parms.Results.samplingRate;
showprogress = parms.Results.showprogress;
roundfreqs = parms.Results.roundfreqs;
saveMatPath = parms.Results.saveMatPath;
MatNameExtraText = parms.Results.MatNameExtraText;
fvector = parms.Results.fvector;
intervals = parms.Results.intervals;
chanID = parms.Results.chanID;
downsampleout = parms.Results.downsampleout;


%Channel restrict
if ~isempty(chanID)
    usechannel = ismember(lfp.channels,chanID);
    lfp.data = lfp.data(:,usechannel);
    lfp.channels = lfp.channels(usechannel);
end

%lfp input
if isstruct(lfp)
    samplingRate = lfp.samplingRate;
elseif isempty(lfp)
    wavespec = lfp;
    return
elseif isnumeric(lfp)
    data_temp = lfp;
    clear lfp
    lfp.data = data_temp;
    lfp.timestamps = [1:length(lfp.data)]'./samplingRate;
end

si = 1./samplingRate;

%Restrict to intervals, with overhang to remove edge effects at transitions
%(then remove later)
overhang = (ncyc)./frange(1);
overint = bsxfun(@(X,Y) X+Y,intervals,overhang.*[-1 1]);
keepIDX = InIntervals(lfp.timestamps,overint);
lfp.data = lfp.data(keepIDX,:);
lfp.timestamps = lfp.timestamps(keepIDX);

%%
if ~isa(lfp.data,'single') || ~isa(lfp.data,'double')
    lfp.data = single(lfp.data);
end

%Frequencies
if ~isempty(fvector)
    freqs = fvector;
else
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
end

if roundfreqs
    freqs = unique(round(freqs));
end

%Filter with wavelets
nfreqs = size(freqs,2);
nchan = size(lfp.data,2);
ntime = ceil(size(lfp.data,1)./downsampleout);
wavespec.data = nan(ntime,nfreqs,nchan);
wavespec.timestamps = downsample(lfp.timestamps,downsampleout);
for cidx = 1:nchan
    for f_i = 1:nfreqs
        if showprogress
            bz_Counter(f_i,nfreqs,'Wavelet Frequency')
        end
        wavelet = MorletWavelet(freqs(f_i),ncyc,si);
         wavespec.data(:,f_i,cidx) = ...
             downsample(FConv(wavelet',lfp.data(:,cidx)),downsampleout);
    end
end

%% Output in buzcode format
%Remove the overhang from intervals
keepIDX = InIntervals(wavespec.timestamps,intervals);
wavespec.data = wavespec.data(keepIDX,:);
wavespec.timestamps = wavespec.timestamps(keepIDX);

wavespec.freqs = freqs;
wavespec.nfreqs = nfreqs;
wavespec.samplingRate = samplingRate./downsampleout;
if isstruct(lfp) && isfield(lfp,'channels')
    wavespec.channels = lfp.channels;
end
wavespec.filterparms.ncyc = ncyc;
wavespec.filterparms.nfreqs = nfreqs;
wavespec.filterparms.frange = frange;
wavespec.filterparms.space = space;

clear lfp

if saveMatPath
    baseName = bz_BasenameFromBasepath(saveMatPath);
    if ~isempty(MatNameExtraText)
        lfpfilename = fullfile(saveMatPath,[baseName,'.wavespec' MatNameExtraText '.lfp.mat']);    
    else
        lfpfilename = fullfile(saveMatPath,[baseName,'.wavespec.lfp.mat']);
    end
    
    s = whos('wavespec');
    if s.bytes>=1073741824%if greater than 2GB
        disp('wavespec variable greater than 2GB, saving as v7.3 .mat file')
        save(lfpfilename,'wavespec','-v7.3')
    else
        save(lfpfilename,'wavespec')
    end
end

end

